/*
 * H4 Chain: General N-Basis RHF, FCI, and VQE at 8 Qubits
 *
 * A linear chain of 4 hydrogen atoms (STO-3G, 1.4 bohr spacing) is a known
 * benchmark for electron correlation and quantum simulation.
 *
 * NOTE: this example's VQE section does not reach the same near-exact
 * convergence eg_42_h2_vqe.c gets for H2. That's expected since8 qubits means
 * 48 ansatz parameters (vs H2's 12), a harder non-convex optimization landscape
 * for simple coordinate-descent optimizer in vqe.c, with real seed sensitivity.
 * This is well-known phenomenon in VQE research (barren plateaus / local minima
 * getting harder at larger qubit count), not something specific to this
 * implementation.
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/second_quant.h"
#include "../physics/vqe.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > H4 Chain: General RHF, FCI, and VQE at 8 Qubits\n\n");

  int n = 4;
  double R = 1.4;
  double centers_pos[4][3];
  basis_function_t *funcs[4];
  double charge[4];

  for (int i = 0; i < n; i++) {
    centers_pos[i][0] = 0;
    centers_pos[i][1] = 0;
    centers_pos[i][2] = i * R;
    funcs[i] = molint_basis_sto3g_h(centers_pos[i]);
    charge[i] = 1.0;
  }
  molecule_t *mol = molecule_alloc(n, charge, centers_pos);

  printf("Geometry: 4 hydrogen atoms, linear chain, %.2f bohr spacing\n\n", R);

  printf(
      "Step 1: General N-basis-function RHF (from physics/molecular_hf.c)\n\n");
  molecular_hf_result_t *hf = molecular_rhf(funcs, n, mol, 4, 1e-10, 200);
  printf("  RHF total energy: %.6f Hartree (converged in %d iterations)\n\n",
         hf->total_energy, hf->iterations);

  printf("Step 2: AO -> MO transform, Jordan-Wigner molecular Hamiltonian "
         "(4 spatial orbitals -> 8 spin-orbitals -> 8 qubits, 256x256)\n\n");

  cmatrix_t *Hcore = molecular_core_hamiltonian(funcs, n, mol);
  double *eri_ao = molecular_eri_tensor(funcs, n);
  double nuclear_repulsion = molecule_nuclear_repulsion(mol);
  double *h_mo = malloc((size_t)n * n * sizeof(double));
  double *eri_mo = malloc((size_t)n * n * n * n * sizeof(double));
  molecular_ao_to_mo(Hcore, eri_ao, hf->C, n, h_mo, eri_mo);
  cmatrix_t *H = second_quant_build_molecular_hamiltonian(n, h_mo, eri_mo,
                                                          nuclear_repulsion);

  printf("Step 3: Exact diagonalization (FCI within this basis)\n\n");
  int dim = H->nrows;
  cmatrix_t *H_copy = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    H_copy->data[i] = H->data[i];
  }

  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  double E_fci = eig->eigenvalues[0];

  double N_expect = 0.0;
  for (int state = 0; state < dim; state++) {
    complex_t amp = CMAT(eig->eigenvectors, state, 0);
    double p = amp.re * amp.re + amp.im * amp.im;
    N_expect += p * __builtin_popcount((unsigned)state);
  }
  printf("  FCI ground state: %.6f Hartree  (<N electrons> = %.4f)\n", E_fci,
         N_expect);
  printf("  Correlation energy: %.6f Hartree\n\n", hf->total_energy - E_fci);

  printf("Step 4: VQE (8 qubits, hardware-efficient ansatz)\n\n");
  vqe_result_t vqe_res = vqe_run(8, 6, H, 30, 0.6, 20260810ULL);
  printf("  VQE converged energy: %.6f Hartree\n\n", vqe_res.energy);

  printf("=== Summary ===\n\n");
  printf("  %-25s %-14.6f\n", "RHF (mean-field)", hf->total_energy);
  printf("  %-25s %-14.6f\n", "FCI (exact, this basis)", E_fci);
  printf("  %-25s %-14.6f\n", "VQE (variational)", vqe_res.energy);

  free(vqe_res.theta_opt);
  eigen_free(eig);
  cmatrix_free(H_copy);
  cmatrix_free(H);
  cmatrix_free(Hcore);
  free(eri_ao);
  free(h_mo);
  free(eri_mo);
  molecular_hf_result_free(hf);
  molecule_free(mol);
  for (int i = 0; i < n; i++) {
    basis_function_free(funcs[i]);
  }

  return 0;
}
