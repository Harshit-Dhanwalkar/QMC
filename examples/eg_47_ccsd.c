/*
 * CCSD (Coupled Cluster Singles and Doubles) on LiH/STO-3G
 *
 * This example runs RHF -> CCSD on LiH, a real multi-orbital,
 * degenerate-orbital (2px/2py) system, and reports RHF energy,
 * CCSD correlation energy, and total CCSD energy.
 */

#include "../core/matrix.h"
#include "../physics/ccsd.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > CCSD: Electron Correlation Energy for LiH/STO-3G\n\n");

  double R = 3.015; // bohr, near-equilibrium LiH bond length
  double c_li[3] = {0.0, 0.0, 0.0};
  double c_h[3] = {0.0, 0.0, R};

  basis_function_t *li_orbs[5];
  molint_basis_sto3g_li(c_li, li_orbs);
  basis_function_t *h_orb = molint_basis_sto3g_h(c_h);

  basis_function_t *basis[6] = {li_orbs[0], li_orbs[1], li_orbs[2],
                                li_orbs[3], li_orbs[4], h_orb};
  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);

  printf("Step 1: RHF (6 basis functions, 4 electrons)\n\n");
  molecular_hf_result_t *hf = molecular_rhf(basis, 6, mol, 4, 1e-12, 200);
  printf("  RHF total energy: %.10f Hartree (%d iterations)\n\n",
         hf->total_energy, hf->iterations);

  printf("Step 2: transform integrals to the MO basis, run CCSD\n\n");
  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 6, mol);
  double *eri_ao = molecular_eri_tensor(basis, 6);
  double *h_mo = malloc(36 * sizeof(double));
  double *eri_mo = malloc(6 * 6 * 6 * 6 * sizeof(double));
  molecular_ao_to_mo(h_ao, eri_ao, hf->C, 6, h_mo, eri_mo);

  ccsd_result_t *ccsd = ccsd_run(6, h_mo, eri_mo, hf->orbital_energies, 4, 0,
                                 hf->total_energy, 1e-10, 100);

  if (!ccsd || !ccsd->converged) {
    fprintf(stderr, "CCSD failed to converge.\n");
  } else {
    printf("  CCSD converged in %d iterations\n", ccsd->iterations);
    printf("  CCSD correlation energy: %.10f Hartree\n",
           ccsd->correlation_energy);
    printf("  CCSD total energy:       %.10f Hartree\n\n", ccsd->total_energy);
    printf("  Correlation recovers %.4f%% additional binding beyond RHF.\n",
           100.0 * fabs(ccsd->correlation_energy) / fabs(hf->total_energy));
    free(ccsd);
  }

  free(h_mo);
  free(eri_mo);
  free(eri_ao);
  cmatrix_free(h_ao);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li_orbs[i]);
  }
  basis_function_free(h_orb);
  molecule_free(mol);
  molecular_hf_result_free(hf);

  return 0;
}
