/*
 * Kohn-Sham LDA DFT for Atoms
 *
 * NOTE: physics/dft.c is a closed-shell Kohn-Sham DFT solver in Local Density
 * Approximation (Slater exchange + Perdew-Zunger 1981 correlation), for same
 * s-orbitals-only scope as hartree_fock.c's RHF same closed-shell
 * doubly-occupied-orbital model. The only physical difference: RHF's exact but
 * expensive nonlocal exchange operator is replaced by a cheap *local*
 * exchange-correlation potential V_xc(n(r)), a function only of density at that
 * point.
 *
 * This example runs He (1 orbital) and Be (2 orbitals: 1s, 2s) through both
 * hartree_fock_atom_s_orbitals and dft_lda_atom_s_orbitals on identical grids,
 * side by side, to make RHF-vs-LDA comparison concrete: expect exact < HF < LDA
 * in total energy magnitude
 */

#include "../core/utils.h"
#include "../physics/dft.h"
#include "../physics/hartree_fock.h"
#include <stdio.h>
#include <stdlib.h>

static void run_atom(const char *name, double Z, int n_orbitals,
                     double exact_ref) {
  int N = 200;
  double r_min = 1e-4, r_max = 14.0;
  double *r = linspace(r_min, r_max, N);

  printf("=== %s (Z=%.0f, %d orbital%s, %d electrons) ===\n\n", name, Z,
         n_orbitals, n_orbitals > 1 ? "s" : "", 2 * n_orbitals);

  hf_result_t *hf =
      hartree_fock_atom_s_orbitals(r, N, Z, n_orbitals, 0.3, 1e-9, 200);
  dft_result_t *dft =
      dft_lda_atom_s_orbitals(r, N, Z, n_orbitals, 0.3, 1e-9, 200);

  printf("  %-20s %-14s %-10s\n", "Method", "E_total (Ha)", "converged");
  printf("  %-20s %-14.6f %-10s\n", "Exact (nonrel.)", exact_ref, "-");
  if (hf) {
    printf("  %-20s %-14.6f %-10s\n", "Hartree-Fock", hf->total_energy,
           hf->converged ? "yes" : "no");
  }

  if (dft) {
    printf("  %-20s %-14.6f %-10s\n", "Kohn-Sham LDA", dft->total_energy,
           dft->converged ? "yes" : "no");
    printf("    (E_hartree = %.6f, E_xc = %.6f)\n", dft->E_hartree, dft->E_xc);
  }
  printf("\n");

  if (hf) {
    hf_result_free(hf);
  }

  if (dft) {
    dft_result_free(dft);
  }

  free(r);
}

int main(void) {
  printf(" > Kohn-Sham LDA DFT for Atoms (Slater exchange + PZ81 "
         "correlation)\n\n");

  run_atom("Helium", 2.0, 1, -2.903724);
  run_atom("Beryllium", 4.0, 2, -14.667356);

  printf("Reference values: exact nonrelativistic (Pekeris/CI-quality) and "
         "RHF-limit numbers from standard references; converged-basis-set LDA "
         "(VWN parametrization) gives He=-2.834836, Be=-14.447209 "
         "(arXiv:2202.00647)\n");

  return 0;
}
