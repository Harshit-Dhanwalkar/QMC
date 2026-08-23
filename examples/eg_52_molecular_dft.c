/*
 * Molecular Multi-Atom Kohn-Sham LDA DFT via Becke Grids: H2 and LiH
 *
 * physics/dft.c implements KS-LDA for single atoms, exploiting spherical
 * symmetry to reduce the problem to a 1D radial grid. This example uses
 * physics/molecular_dft.c's generalization to arbitrary molecules : a 3D Becke
 * fuzzy-cell atomic-partition grid, with no symmetry assumption. It runs H2 and
 * LiH, and for each also runs plain RHF for comparison, showing how (and how
 * much) LDA's exchange- correlation treatment differs from Hartree-Fock's
 * exact-exchange-only treatment for these small systems.
 */

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
  molecular_dft_result_t *dft =
      molecular_ks_lda_default(basis, n_basis, mol, n_electrons, grid);

  if (!dft || !dft->converged) {
    fprintf(stderr, "  KS-LDA failed to converge.\n");
  } else {
    printf("  KS-LDA total energy:  %.8f Hartree (%d iterations)\n",
           dft->total_energy, dft->iterations);
    printf("  Energy breakdown: core=%.6f  coulomb=%.6f  xc=%.6f  "
           "nuclear=%.6f\n",
           dft->e_core, dft->e_coulomb, dft->e_xc, dft->e_nuclear);
    printf("  RHF - KS-LDA = %.6f Hartree (LDA's local-density XC treatment "
           "vs. RHF's exact nonlocal exchange : neither is strictly 'more "
           "correct' for a minimal STO-3G basis)\n\n",
           hf->total_energy - dft->total_energy);
  }

  if (dft) {
    molecular_dft_result_free(dft);
  }
  molecular_grid_free(grid);
  molecular_hf_result_free(hf);
}

int main(void) {
  printf(" > Molecular Kohn-Sham LDA DFT (Becke Grid): H2 and LiH\n\n");

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
