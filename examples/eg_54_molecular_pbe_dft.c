/*
 * Molecular Multi-Atom Kohn-Sham PBE (GGA) DFT via Becke Grids: H2 and LiH
 *
 * Run KS-LDA (Refernce: Perdew-Zunger 1981) on Becke fuzzy-cell 3D grid.
 * This example runs physics/dft.c's PBE96 GGA functional (Refernces: Perdew,
 * Burke & Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996)) instead, via
 * physics/molecular_dft.c's molecular_ks_pbe. PBE adds a gradient-dependent
 * correction on top of LDA (built on Perdew-Wang 1992 correlation rather than
 * Perdew-Zunger, since that is what PBE's enhancement factor is analytically
 * parametrized against), so it needs each basis function's gradient in addition
 * to its value at every grid point : that's molecular_integrals.h's
 * basis_function_gradient.
 * Runs RHF, KS-LDA and KS-PBE side by side for H2 and LiH so the three levels
 * of exchange-correlation treatment can be compared directly on the same
 * minimal STO-3G basis.
 */

#include "../physics/dft.h"
#include "../physics/molecular_dft.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <stdio.h>

static void run_comparison(const char *label, basis_function_t **basis,
                           int n_basis, const molecule_t *mol,
                           int n_electrons) {
  printf(" > %s\n\n", label);

  molecular_hf_result_t *hf =
      molecular_rhf(basis, n_basis, mol, n_electrons, 1e-10, 200);
  printf("  RHF total energy:     %.8f Hartree (%d iterations)\n",
         hf->total_energy, hf->iterations);

  molecular_grid_t *grid = molecular_grid_build_default(mol);

  molecular_dft_result_t *lda =
      molecular_ks_lda_default(basis, n_basis, mol, n_electrons, grid);
  molecular_dft_result_t *pbe =
      molecular_ks_pbe_default(basis, n_basis, mol, n_electrons, grid);

  if (!lda || !lda->converged) {
    fprintf(stderr, "  KS-LDA failed to converge.\n");
  } else {
    printf("  KS-LDA total energy:  %.8f Hartree (%d iterations)\n",
           lda->total_energy, lda->iterations);
  }

  if (!pbe || !pbe->converged) {
    fprintf(stderr, "  KS-PBE failed to converge.\n");
  } else {
    printf("  KS-PBE total energy:  %.8f Hartree (%d iterations)\n",
           pbe->total_energy, pbe->iterations);
    printf("  Energy breakdown: core=%.6f  coulomb=%.6f  xc=%.6f  "
           "nuclear=%.6f\n",
           pbe->e_core, pbe->e_coulomb, pbe->e_xc, pbe->e_nuclear);
  }

  if (lda && lda->converged && pbe && pbe->converged) {
    printf("  KS-LDA - KS-PBE = %.6f Hartree (PBE's gradient correction on "
           "top of LDA's local-density treatment)\n",
           lda->total_energy - pbe->total_energy);
    printf("  RHF - KS-PBE    = %.6f Hartree\n\n",
           hf->total_energy - pbe->total_energy);
  }

  if (lda) {
    molecular_dft_result_free(lda);
  }
  if (pbe) {
    molecular_dft_result_free(pbe);
  }
  molecular_grid_free(grid);
  molecular_hf_result_free(hf);
}

int main(void) {
  printf(" > Molecular Kohn-Sham PBE (GGA) DFT (Becke Grid): H2 and LiH\n\n");

  {
    double R = 1.4;
    double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
    basis_function_t *h0 = molint_basis_sto3g_h(c0);
    basis_function_t *h1 = molint_basis_sto3g_h(c1);
    basis_function_t *basis[2] = {h0, h1};
    const double charges[2] = {1.0, 1.0};
    double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
    molecule_t *mol = molecule_alloc(2, charges, centers);

    run_comparison("H2/STO-3G, R=1.4 bohr", basis, 2, mol, 2);

    basis_function_free(h0);
    basis_function_free(h1);
    molecule_free(mol);
  }

  {
    double R = 3.015;
    double c_li[3] = {0, 0, 0}, c_h[3] = {0, 0, R};
    basis_function_t *li_orbs[5];
    molint_basis_sto3g_li(c_li, li_orbs);
    basis_function_t *h_orb = molint_basis_sto3g_h(c_h);
    basis_function_t *basis[6] = {li_orbs[0], li_orbs[1], li_orbs[2],
                                  li_orbs[3], li_orbs[4], h_orb};
    const double charges[2] = {3.0, 1.0};
    double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
    molecule_t *mol = molecule_alloc(2, charges, centers);

    run_comparison("LiH/STO-3G, R=3.015 bohr", basis, 6, mol, 4);

    for (int i = 0; i < 5; i++) {
      basis_function_free(li_orbs[i]);
    }
    basis_function_free(h_orb);
    molecule_free(mol);
  }

  return 0;
}
