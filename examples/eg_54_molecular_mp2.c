/*
 * General Molecular MP2: Beyond the s-Orbital Restriction
 *
 * NOTE: eg_32_mp2.c's MP2 is restricted to s-only atomic orbitals, so it can
 * only ever recover the l=0 slice of the true correlation energy - there's no
 * way to describe correlation into p, d, ... virtuals with that machinery.
 * molecular_mp2() : it works on the same general Gaussian-basis MO integrals
 * ccsd_run() uses (McMurchie-Davidson AO integrals + molecular_ao_to_mo), so
 * real molecules with real angular momentum in their virtual space are fully in
 * scope.
 *
 * This runs RHF -> MP2 on LiH/STO-3G (which has genuine 2px/2py virtual
 * character on Li), both with and without a frozen Li 1s core.
 */

#include "../core/matrix.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/mp2.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > General Molecular MP2: Correlation Energy for LiH/STO-3G\n\n");

  double R = 3.015; // bohr, same near-equilibrium bond length as eg_47_ccsd
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

  printf("Step 2: transform integrals to the MO basis, run MP2\n\n");
  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 6, mol);
  double *eri_ao = molecular_eri_tensor(basis, 6);
  double *h_mo = malloc(36 * sizeof(double));
  double *eri_mo = malloc(6 * 6 * 6 * 6 * sizeof(double));
  molecular_ao_to_mo(h_ao, eri_ao, hf->C, 6, h_mo, eri_mo);

  molecular_mp2_result_t full =
      molecular_mp2(6, eri_mo, hf->orbital_energies, 4, 0, hf->total_energy);
  printf("  All-electron MP2:\n");
  printf("    correlation energy: %.10f Hartree\n", full.e_mp2);
  printf("    total energy:       %.10f Hartree\n", full.e_total);
  printf("    (reference: -0.0128683238 Hartree, |diff| = %.2e)\n\n",
         fabs(full.e_mp2 - (-0.012868323831743156)));

  molecular_mp2_result_t frozen =
      molecular_mp2(6, eri_mo, hf->orbital_energies, 4, 1, hf->total_energy);
  printf("  Frozen-core MP2 (Li 1s excluded):\n");
  printf("    correlation energy: %.10f Hartree\n", frozen.e_mp2);
  printf("    total energy:       %.10f Hartree\n", frozen.e_total);
  printf("    (reference: -0.0126403353 Hartree, |diff| = %.2e)\n\n",
         fabs(frozen.e_mp2 - (-0.012640335346251669)));

  printf("  Freezing the Li 1s core changes the correlation energy by only "
         "%.2e Hartree (%.2f%%): 1s orbital is too tightly bound and too far "
         "in energy from the virtuals to contribute much beyond-HF "
         "correlation, exactly as chemical intuition expects.\n",
         fabs(full.e_mp2 - frozen.e_mp2),
         100.0 * fabs(full.e_mp2 - frozen.e_mp2) / fabs(full.e_mp2));

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
