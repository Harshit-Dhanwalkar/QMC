/*
 * CCSD(T) (Perturbative Triples Correction) on an Asymmetric H4 Cluster
 */

#include "../core/matrix.h"
#include "../physics/ccsd_t.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > CCSD(T): Perturbative Triples on an Asymmetric H4 Cluster\n\n");

  // Arbitrary non-collinear geometry (bohr)
  double centers[4][3] = {
      {0.0, 0.0, 0.0}, {0.9, 0.3, 0.1}, {1.7, -0.4, 0.6}, {2.9, 0.5, -0.3}};
  const double charges[4] = {1.0, 1.0, 1.0, 1.0};

  basis_function_t *basis[4];
  for (int i = 0; i < 4; i++) {
    basis[i] = molint_basis_sto3g_h(centers[i]);
  }

  molecule_t *mol = molecule_alloc(4, charges, centers);

  printf("Step 1: RHF (4 basis functions, 4 electrons)\n\n");
  molecular_hf_result_t *hf = molecular_rhf(basis, 4, mol, 4, 1e-12, 200);
  printf("  RHF total energy:      %.10f Hartree\n\n", hf->total_energy);

  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 4, mol);
  double *eri_ao = molecular_eri_tensor(basis, 4);
  double *h_mo = malloc(16 * sizeof(double));
  double *eri_mo = malloc(4 * 4 * 4 * 4 * sizeof(double));

  molecular_ao_to_mo(h_ao, eri_ao, hf->C, 4, h_mo, eri_mo);

  printf("Step 2: CCSD(T) (CCSD to convergence, then the perturbative "
         "triples correction)\n\n");
  ccsdt_result_t *res = ccsdt_run(4, h_mo, eri_mo, hf->orbital_energies, 4, 0,
                                  hf->total_energy, 1e-10, 100);

  if (!res) {
    fprintf(stderr, "CCSD(T) failed to converge.\n");
  } else {
    printf("  CCSD correlation energy:      %.10f Hartree (%d iterations)\n",
           res->ccsd_correlation_energy, res->ccsd_iterations);
    printf("  (T) perturbative correction:  %.10f Hartree\n",
           res->perturbative_correction);
    printf("  CCSD(T) total energy:         %.10f Hartree\n\n",
           res->total_energy);
    printf("  (T) correction is %.3f%% of the CCSD correlation energy on this "
           "non-symmetric cluster : small but nonzero, as expected for a "
           "perturbative refinement.\n",
           100.0 * res->perturbative_correction / res->ccsd_correlation_energy);
    free(res);
  }

  free(h_mo);
  free(eri_mo);
  free(eri_ao);
  cmatrix_free(h_ao);
  for (int i = 0; i < 4; i++) {
    basis_function_free(basis[i]);
  }
  molecule_free(mol);
  molecular_hf_result_free(hf);

  return 0;
}
