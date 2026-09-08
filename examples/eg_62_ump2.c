/*
 * UMP2: Post-HF Correlation Energy for Open-Shell Radicals
 *
 * NOTE: eg_59_uhf_openshell.c gets UHF energies for radicals, but stops at the
 * mean-field level - there is no way to add a correlation correction for an
 * open-shell system the way molecular_mp2() does for RHF's closed-shell
 * reference.
 *
 * Runs UHF -> UMP2 on the lithium atom's doublet ground state (2 \alpha, 1
 * \beta), the same system eg_59 uses, and checks against an independent
 * reference.
 */

#include "../core/matrix.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/mp2.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > UMP2: Correlation Energy for the Lithium Doublet Radical\n\n");

  double c_li[3] = {0.0, 0.0, 0.0};
  basis_function_t *li_orbs[5];
  molint_basis_sto3g_li(c_li, li_orbs);

  const double charges[1] = {3.0};
  double centers[1][3] = {{0.0, 0.0, 0.0}};
  molecule_t *mol = molecule_alloc(1, charges, centers);

  printf("  Step 1: UHF (5 basis functions, 2 \\alpha + 1 \\beta electrons)\n\n");
  molecular_uhf_result_t *uhf =
      molecular_uhf(li_orbs, 5, mol, 2, 1, 1e-12, 300);
  printf("    UHF total energy: %.10f Hartree (%d iterations)\n",
         uhf->total_energy, uhf->iterations);
  printf("    <S^2> = %.6f (exact 0.75 for an uncontaminated doublet)\n\n",
         uhf->spin_squared);

  printf("  Step 2: transform integrals to the alpha and beta MO bases, run UMP2\n\n");
  cmatrix_t *h_ao = molecular_core_hamiltonian(li_orbs, 5, mol);
  double *eri_ao = molecular_eri_tensor(li_orbs, 5);

  double *eri_aaaa = malloc(5 * 5 * 5 * 5 * sizeof(double));
  double *eri_bbbb = malloc(5 * 5 * 5 * 5 * sizeof(double));
  double *eri_aabb = malloc(5 * 5 * 5 * 5 * sizeof(double));
  double *h_mo_scratch = malloc(25 * sizeof(double));

  molecular_ao_to_mo(h_ao, eri_ao, uhf->C_alpha, 5, h_mo_scratch, eri_aaaa);
  molecular_ao_to_mo(h_ao, eri_ao, uhf->C_beta, 5, h_mo_scratch, eri_bbbb);
  molecular_ao_to_mo_eri_mixed(eri_ao, uhf->C_alpha, uhf->C_beta, 5, eri_aabb);

  molecular_ump2_result_t ump2 = molecular_ump2(
      5, eri_aaaa, eri_bbbb, eri_aabb, uhf->orbital_energies_alpha,
      uhf->orbital_energies_beta, 2, 1, 0, 0, uhf->total_energy);

  printf("    Same-spin (alpha-alpha):     %.10f Hartree\n", ump2.e_aa);
  printf("    Same-spin (beta-beta):       %.10f Hartree\n", ump2.e_bb);
  printf("    Opposite-spin (alpha-beta):  %.10f Hartree\n", ump2.e_ab);
  printf("    Total UMP2 correlation:      %.10f Hartree\n", ump2.e_mp2);
  printf("    UHF+UMP2 total energy:       %.10f Hartree\n\n", ump2.e_total);

  printf("    (Independent reference: UHF=-7.3155259813, UMP2 corr=-0.0002564094 "
         "Hartree)\n");
  printf("    |UHF diff| = %.2e, |UMP2 corr diff| = %.2e\n\n",
         fabs(uhf->total_energy - (-7.3155259813)),
         fabs(ump2.e_mp2 - (-0.0002564094444800904)));

  free(h_mo_scratch);
  free(eri_aabb);
  free(eri_bbbb);
  free(eri_aaaa);
  free(eri_ao);
  cmatrix_free(h_ao);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li_orbs[i]);
  }
  molecule_free(mol);
  molecular_uhf_result_free(uhf);

  return 0;
}
