/*
 * Unitary Coupled Cluster (UCC) Ansatz, CCSD-Seeded, on H2/STO-3G
 *
 * UCC is implemented and validated in tests/test_ucc.c, but had no dedicated
 * example. Unlike CCSD's plain \exp(T) (not unitary), UCC's \exp(\kappa) with
 * anti-Hermitian \kappa is unitary by construction : a legitimate ansatz for
 * real quantum circuit. This example wires the full pipeline end to end:
 *
 *   RHF -> MO integrals -> CCSD (for amplitudes) -> UCC generator ->
 *   state preparation -> energy expectation vs. exact FCI
 *
 * and shows that seeding UCC's parameters directly from CCSD's converged
 * amplitudes already lands (for this 2-electron system) essentially exactly on
 * FCI ground state : physical basis for "UCCSD" as a practical VQE ansatz
 * initialization strategy.
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/ccsd.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/second_quant.h"
#include "../physics/ucc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > UCC Ansatz: CCSD-Seeded State Preparation for H2/STO-3G\n\n");

  double R = 1.4;
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};
  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);

  printf("Step 1: RHF + MO integrals + CCSD (for amplitudes)\n\n");
  molecular_hf_result_t *hf = molecular_rhf(basis, 2, mol, 2, 1e-12, 200);
  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 2, mol);
  double *eri_ao = molecular_eri_tensor(basis, 2);
  double *h_mo = malloc(4 * sizeof(double));
  double *eri_mo = malloc(16 * sizeof(double));
  molecular_ao_to_mo(h_ao, eri_ao, hf->C, 2, h_mo, eri_mo);

  ccsd_amplitudes_t *amp = NULL;
  ccsd_result_t *ccsd = ccsd_run_ex(2, h_mo, eri_mo, hf->orbital_energies, 2, 0,
                                    hf->total_energy, 1e-12, 100, &amp);
  printf("  CCSD total energy: %.10f Hartree\n\n", ccsd->total_energy);

  printf("Step 2: build the UCC generator from CCSD's converged T1/T2 "
         "amplitudes\n\n");
  ucc_single_t *singles = NULL;
  double *theta_s = NULL;
  int n_singles = 0;
  ucc_double_t *doubles = NULL;
  double *theta_d = NULL;
  int n_doubles = 0;
  ucc_excitations_from_ccsd_amplitudes(amp, &singles, &theta_s, &n_singles,
                                       &doubles, &theta_d, &n_doubles);
  printf("  %d singles, %d doubles excitation(s) extracted from CCSD "
         "amplitudes\n",
         n_singles, n_doubles);

  int n_modes = amp->nso; // 4 spin orbitals for H2/STO-3G
  cmatrix_t *generator = ucc_build_generator(
      n_modes, singles, theta_s, n_singles, doubles, theta_d, n_doubles);

  /* Hartree-Fock reference determinant: modes 0,1 occupied (two lowest spin
   * orbitals), big-endian bit convention matching second_quant.c. */
  cvector_t *hf_ref = cvector_alloc(1 << n_modes);
  for (int k = 0; k < (1 << n_modes); k++) {
    hf_ref->data[k] = c_zero();
  }

  hf_ref->data[0b1100] = c_real(1.0);

  cvector_t *psi = ucc_prepare_state(generator, hf_ref);

  printf("\nStep 3: energy expectation "
         "<\\psi(\\theta_{CCSD})|H|\\psi(\\theta_{CCSD})> vs. exact FCI\n\n");
  cmatrix_t *H = second_quant_build_molecular_hamiltonian(
      2, h_mo, eri_mo, molecule_nuclear_repulsion(mol));

  double e_ucc = 0.0;

  int dim = 1 << n_modes;
  for (int r = 0; r < dim; r++) {
    complex_t Hpsi_r = c_zero();

    for (int c = 0; c < dim; c++) {
      Hpsi_r = c_add(Hpsi_r, c_mul(CMAT(H, r, c), psi->data[c]));
    }

    e_ucc += c_mul(c_conj(psi->data[r]), Hpsi_r).re;
  }

  eigen_t *eig = cmatrix_eigh_complex(H);
  double e_fci = eig->eigenvalues[0];
  for (int i = 1; i < eig->n; i++) {
    if (eig->eigenvalues[i] < e_fci) {
      e_fci = eig->eigenvalues[i];
    }
  }

  printf("  UCC (CCSD-seeded) energy: %.10f Hartree\n", e_ucc);
  printf("  Exact FCI ground energy:  %.10f Hartree\n", e_fci);
  printf("  Difference:                %.2e Hartree\n\n", fabs(e_ucc - e_fci));
  printf("  UCC is variational (energy >= FCI always); CCSD-seeded parameters "
         "already essentially saturate that bound for this 2-electron system, "
         "justification for using CCSD amplitudes as a VQE UCCSD "
         "initialization.\n");

  cvector_free(psi);
  cvector_free(hf_ref);
  cmatrix_free(generator);
  cmatrix_free(H);
  eigen_free(eig);
  free(singles);
  free(theta_s);
  free(doubles);
  free(theta_d);
  ccsd_amplitudes_free(amp);
  free(ccsd);
  free(h_mo);
  free(eri_mo);
  free(eri_ao);
  cmatrix_free(h_ao);
  basis_function_free(h0);
  basis_function_free(h1);
  molecule_free(mol);
  molecular_hf_result_free(hf);

  return 0;
}
