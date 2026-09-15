/*
 * Linear-Response TDDFT (TDA): Electronic Excitation Energies from a
 * Ground-State DFT Calculation
 *
 * NOTE: molecular_ks_lda (physics/molecular_dft.c) gives only the ground state.
 * TDDFT reuses that converged calculation's orbitals/energies to additionally
 * compute vertical excitation energies - energies of light absorbed in
 * promoting an electron from an occupied to a virtual Kohn-Sham orbital,
 * corrected for electron-electron interaction via same exchange-correlation
 * functional. This is the standard first-principles route to UV-Vis absorption
 * spectra
 */

#include "../physics/molecular_dft.h"
#include "../physics/molecular_integrals.h"
#include "../physics/tddft.h"
#include <stdio.h>

static void run_tddft(const char *label, basis_function_t **basis, int n_basis,
                      const molecule_t *mol, int n_electrons) {
  printf(" > %s\n\n", label);

  molecular_grid_t *grid = molecular_grid_build_default(mol);
  molecular_dft_result_t *ks =
      molecular_ks_lda_default(basis, n_basis, mol, n_electrons, grid);

  printf("  Ground state KS-LDA(PZ81) energy: %.8f Hartree\n",
         ks->total_energy);

  tddft_result_t *td = tddft_tda_lda(basis, n_basis, mol, ks, grid);

  printf("  Active space: %d occupied x %d virtual orbitals (%d possible "
         "single excitations)\n\n",
         td->n_occ, td->n_virt, td->dim);

  printf("  Vertical singlet excitation energies (TDA-LDA):\n");
  for (int k = 0; k < td->dim; k++) {
    double ev = td->excitation_energies[k] * 27.211386245988;

    printf("    S%d: %.6f Hartree  (%.4f eV)\n", k + 1,
           td->excitation_energies[k], ev);
  }
  printf("\n");

  tddft_result_free(td);
  molecular_dft_result_free(ks);
  molecular_grid_free(grid);
}

int main(void) {
  printf(" > TDDFT: Electronic Excitation Energies\n\n");

  {
    double R = 1.4;
    double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
    basis_function_t *h0 = molint_basis_sto3g_h(c0);
    basis_function_t *h1 = molint_basis_sto3g_h(c1);
    basis_function_t *basis[2] = {h0, h1};
    const double charges[2] = {1.0, 1.0};
    double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
    molecule_t *mol = molecule_alloc(2, charges, centers);

    run_tddft("H2/STO-3G, R=1.4 bohr", basis, 2, mol, 2);

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

    run_tddft("LiH/STO-3G, R=3.015 bohr", basis, 6, mol, 4);

    for (int i = 0; i < 5; i++) {
      basis_function_free(li_orbs[i]);
    }
    basis_function_free(h_orb);
    molecule_free(mol);
  }

  printf(" NOTE: LiH's second and third excitations are degenerate (equal "
         "energy) - expected, since they correspond to excitations into a pair "
         "of degenerate pi-symmetry virtual orbitals for a linear diatomic.\n");

  return 0;
}
