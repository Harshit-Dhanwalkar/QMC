/*
 * Unrestricted Hartree-Fock (UHF): Open-Shell Radicals
 *
 * molecular_rhf() forces \alpha and \beta electrons into same spatial orbitals,
 * paired up two at a time. molecular_uhf() \alpha and \beta electrons each get
 * their own density matrix and their own Fock operator, self-consistently
 * coupled only through the shared Coulomb (J) term built from the total density
 * (Reference: UHF equations, Szabo & Ostlund Ch. 3.8)
 *
 * This example runs three cases that build on each other:
 *   1. The hydrogen atom (1 electron - the simplest possible open shell,
 *      no beta electrons at all).
 *   2. The lithium atom's doublet ground state (2 alpha, 1 beta - a real
 *      radical), checked against an independent PySCF reference.
 *   3. Li+ (2 electrons, closed-shell) run through the *same* UHF code
 *      with n_alpha = n_beta, to show that UHF reduces to RHF exactly
 *      when there's nothing asymmetric left to describe.
 */

#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void run_h_atom(void) {
  printf(" === Hydrogen atom (1 electron: 1 \\alpha, 0 \\beta) ===\n\n");

  double center[3] = {0.0, 0.0, 0.0};
  basis_function_t *h_orb = molint_basis_sto3g_h(center);

  const double charge[1] = {1.0};
  double centers[1][3] = {{0.0, 0.0, 0.0}};
  molecule_t *mol = molecule_alloc(1, charge, centers);

  molecular_uhf_result_t *uhf = molecular_uhf(&h_orb, 1, mol, 1, 0, 1e-10, 100);
  if (!uhf) {
    fprintf(stderr, "  UHF failed to allocate for H atom.\n\n");

    molecule_free(mol);
    basis_function_free(h_orb);

    return;
  }

  printf("  converged=%s in %d iterations\n", uhf->converged ? "yes" : "no",
         uhf->iterations);
  printf("  E_total    = %.10f Hartree\n", uhf->total_energy);
  printf("  alpha HOMO = %.10f Hartree (only occupied orbital)\n",
         uhf->orbital_energies_alpha[0]);
  printf("  <S^2>      = %.6f (exact Sz(Sz+1) for a lone electron: 0.75)\n\n",
         uhf->spin_squared);

  molecular_uhf_result_free(uhf);
  molecule_free(mol);
  basis_function_free(h_orb);
}

static molecular_uhf_result_t *run_li(int n_alpha, int n_beta,
                                      basis_function_t *li[5],
                                      const molecule_t *mol) {
  return molecular_uhf(li, 5, mol, n_alpha, n_beta, 1e-10, 300);
}

int main(void) {
  printf(" > Unrestricted Hartree-Fock: Open-Shell Radicals\n\n");

  run_h_atom();

  printf("  === Lithium atom (doublet radical: 2 alpha, 1 beta) ===\n\n");

  double li_center[3] = {0.0, 0.0, 0.0};
  basis_function_t *li[5];
  if (!molint_basis_sto3g_li(li_center, li)) {
    fprintf(stderr, "  Failed to build Li STO-3G basis.\n");

    return 1;
  }

  const double li_charge[1] = {3.0};
  double li_pos[1][3] = {{0.0, 0.0, 0.0}};
  molecule_t *mol = molecule_alloc(1, li_charge, li_pos);

  molecular_uhf_result_t *doublet = run_li(2, 1, li, mol);
  if (!doublet) {
    fprintf(stderr, "  UHF failed to allocate for Li atom.\n");
    molecule_free(mol);

    for (int i = 0; i < 5; i++) {
      basis_function_free(li[i]);
    }

    return 1;
  }

  printf("  converged=%s in %d iterations\n", doublet->converged ? "yes" : "no",
         doublet->iterations);
  printf("  E_total = %.10f Hartree\n", doublet->total_energy);
  printf(
      "  (independent PySCF reference: -7.3155259813 Hartree, |diff| = %.2e)\n",
      fabs(doublet->total_energy - (-7.3155259813)));

  printf("  alpha orbital energies:");
  for (int k = 0; k < doublet->n_basis; k++) {
    printf(" %+.6f", doublet->orbital_energies_alpha[k]);
  }

  printf("\n  beta  orbital energies:");
  for (int k = 0; k < doublet->n_basis; k++) {
    printf(" %+.6f", doublet->orbital_energies_beta[k]);
  }
  printf("\n  <S^2> = %.6f (exact Sz(Sz+1) for a clean doublet: 0.75 - 2s "
         "electron is well-separated in energy from other configurations here, "
         "so there's no spin contamination)\n\n",
         doublet->spin_squared);

  printf("  === Li+ (2 electrons, closed-shell: n_alpha = n_beta = 1) ===\n");

  molecular_uhf_result_t *liplus = run_li(1, 1, li, mol);
  molecular_hf_result_t *rhf_liplus = molecular_rhf(li, 5, mol, 2, 1e-10, 300);

  if (liplus && rhf_liplus) {
    printf("  UHF(Li+) E_total = %.10f Hartree\n", liplus->total_energy);
    printf("  RHF(Li+) E_total = %.10f Hartree (|diff| = %.2e)\n",
           rhf_liplus->total_energy,
           fabs(liplus->total_energy - rhf_liplus->total_energy));
    printf("  <S^2> = %.6f (exact 0 for a closed-shell singlet)\n\n",
           liplus->spin_squared);
  } else {
    fprintf(stderr, "  UHF/RHF failed to allocate for Li+.\n\n");
  }

  molecular_uhf_result_free(liplus);
  molecular_hf_result_free(rhf_liplus);
  molecular_uhf_result_free(doublet);
  molecule_free(mol);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }

  return 0;
}
