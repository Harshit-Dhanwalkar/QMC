/*
 * LiH: p-Orbitals, Frozen-Core Reduction, and VQE at 10 Qubits
 *
 * LiH/STO-3G have p-orbitals (molint_basis_sto3g_li: 1s core + 2s + full 2p
 * shell) and a heteronuclear pair. LiH's full minimal basis is 6 spatial
 * orbitals (12 qubits, 4096-dim), but freezing Li's tightly-bound, chemically
 * inert 1s core, reduces this to 5 active spatial orbitals (10 qubits,
 * 1024-dim) at a cost of only about 0.23 mHartree versus the exact full-space
 * answer.
 */

#include "../core/matrix.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/second_quant.h"
#include "../physics/vqe.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_LAPACK
#include "../core/linalg/complex_eigh.h"
#endif

int main(void) {
  printf(" > LiH: p-Orbitals, Frozen-Core Reduction, and VQE at 10 Qubits\n\n");

  double R = 3.015; // bohr, ~1.596 Angstrom, near LiH's experimental
                    // equilibrium bond length
  double li_center[3] = {0, 0, 0};
  double h_center[3] = {0, 0, R};

  printf("Step 1: STO-3G basis (Li: 1s + 2s + 2p shell, H: 1s)\n\n");
  basis_function_t *li[5];
  molint_basis_sto3g_li(li_center, li);
  basis_function_t *h = molint_basis_sto3g_h(h_center);
  basis_function_t *basis[6] = {li[0], li[1], li[2], li[3], li[4], h};

  double charge[2] = {3.0, 1.0};
  double centers_pos[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charge, centers_pos);

  printf("Step 2: General N-basis RHF (6 basis functions, 4 electrons)\n\n");
  molecular_hf_result_t *hf = molecular_rhf(basis, 6, mol, 4, 1e-9, 300);
  printf("  LiH RHF total energy: %.6f Hartree (converged in %d iterations)\n"
         "  (Theoretical reference: GAMESS/Molpro give -7.8633823379 Hartree "
         "at a similar geometry)\n\n",
         hf->total_energy, hf->iterations);

  printf("Step 3: AO -> MO transform, then freeze Li's 1s core "
         "(6 spatial orbitals -> 5 active, 12 qubits -> 10)\n\n");

  cmatrix_t *Hcore = molecular_core_hamiltonian(basis, 6, mol);
  double *eri_ao = molecular_eri_tensor(basis, 6);
  double nuclear_repulsion = molecule_nuclear_repulsion(mol);
  double *h_mo = malloc(36 * sizeof(double));
  double *eri_mo = malloc(1296 * sizeof(double));
  molecular_ao_to_mo(Hcore, eri_ao, hf->C, 6, h_mo, eri_mo);

  cmatrix_t *H = second_quant_build_molecular_hamiltonian_frozen_core(
      6, 1, h_mo, eri_mo, nuclear_repulsion);
  printf("  Frozen-core Hamiltonian: %dx%d\n\n", H->nrows, H->ncols);

#ifdef USE_LAPACK
  printf("Step 4: Exact diagonalization (FCI within the frozen-core active "
         "space) : USE_LAPACK build, so this is tractable\n\n");
  int dim = H->nrows;
  cmatrix_t *H_copy = cmatrix_alloc(dim, dim);

  for (int i = 0; i < dim * dim; i++) {
    H_copy->data[i] = H->data[i];
  }

  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  double E_fci = eig->eigenvalues[0];
  printf("  Frozen-core FCI: %.6f Hartree  (full 12-qubit FCI, cross-validated "
         ": -7.876960 Hartree-frozen-core error ~%.4f mHartree)\n\n",
         E_fci, 1000.0 * fabs(E_fci - (-7.876960)));

  printf("Step 5: VQE (10 qubits, hardware-efficient ansatz)\n\n");
  vqe_result_t vqe_res = vqe_run(10, 4, H, 15, 0.6, 20260810ULL);
  printf("  VQE converged energy: %.6f Hartree\n\n", vqe_res.energy);

  printf("=== Summary ===\n\n");
  printf("  %-30s %-14.6f\n", "RHF (mean-field)", hf->total_energy);
  printf("  %-30s %-14.6f\n", "Frozen-core FCI (exact)", E_fci);
  printf("  %-30s %-14.6f\n", "VQE (variational)", vqe_res.energy);
  printf("\n  10-qubit VQE run does not reach H2's near-exact convergence with "
         "a optimizer budget : 40 ansatz parameters (10 qubits x 4 layers) is "
         "a harder non-convex landscape than H2's 12, and 10-qubit state "
         "simulation itself is dominant cost here. The Hamiltonian and FCI "
         "energy above are exact regardless.\n");

  free(vqe_res.theta_opt);
  eigen_free(eig);
  cmatrix_free(H_copy);
#else
  printf("Step 4: (skipped : exact diagonalization of a 1024x1024 complex "
         "Hermitian matrix needs USE_LAPACK=1; rebuild with `make "
         "USE_LAPACK=1` to see the FCI energy and VQE run here. RHF and the "
         "Hamiltonian construction above are unaffected either way.)\n\n");
#endif

  cmatrix_free(H);
  cmatrix_free(Hcore);
  free(eri_ao);
  free(h_mo);
  free(eri_mo);
  molecular_hf_result_free(hf);
  molecule_free(mol);
  basis_function_free(h);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }

  return 0;
}
