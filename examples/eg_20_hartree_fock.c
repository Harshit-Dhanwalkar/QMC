/*
 * Restricted Hartree-Fock - Beyond the Effective-Charge Trick
 *
 * SCF loop instead solves self-consistently for orbital's actual radial *shape*
 * - direct + exchange mean-field potentials built from current orbitals, a
 * dense Fock matrix diagonalized each iteration, repeat until orbitals stop
 * changing. Because HF searches a strictly larger space of trial wavefunctions
 * than single-parameter exponential ansatz, it must land at an equal-or-lower
 * (more -ve) energy
 */

#include "../core/utils.h"
#include "../physics/hartree_fock.h"
#include "../physics/helium.h"
#include <stdio.h>
#include <stdlib.h>

static void run_atom(const char *name, double Z, int n_orbitals, double r_max) {
  int N = 160;
  double *r = linspace(1e-4, r_max, N);

  hf_result_t *res =
      hartree_fock_atom_s_orbitals(r, N, Z, n_orbitals, 0.4, 1e-8, 300);

  if (!res) {
    printf("   %s: SCF setup failed\n\n", name);
    free(r);

    return;
  }

  printf("   %s (Z=%.0f, %d doubly-occupied s-orbital%s):\n", name, Z,
         n_orbitals, n_orbitals == 1 ? "" : "s");
  printf("     converged=%s in %d iterations\n", res->converged ? "yes" : "no",
         res->iterations);
  for (int k = 0; k < n_orbitals; k++) {
    printf("     orbital %d energy: %12.6f Hartree\n", k,
           res->orbital_energies[k]);
  }
  printf("     total electronic energy: %12.6f Hartree\n", res->total_energy);

  hf_result_free(res);
  free(r);
  printf("\n");
}

int main(void) {
  printf(" > Restricted Hartree-Fock: s-Orbital Closed-Shell Atoms\n\n");

  double z_eff_he = helium_optimal_zeff_analytic(2.0);
  double e_simple_he = helium_ground_state_energy_analytic(2.0);
  printf("   helium.c's simple product-orbital variational estimate:\n");
  printf("     Z_eff = %.6f, E = %.6f Hartree\n\n", z_eff_he, e_simple_he);

  run_atom("Helium", 2.0, 1, 12.0);
  run_atom("Beryllium", 4.0, 2, 10.0);

  printf("   Helium's HF energy above should be <= %.6f (estimate)\n",
         e_simple_he);

  return 0;
}
