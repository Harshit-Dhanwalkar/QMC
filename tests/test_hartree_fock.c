/*
 * Test + demonstration: restricted Hartree-Fock (closed-shell, s-orbitals
 * only), via physics/hartree_fock.c.
 *
 * 1. Helium (Z=2, n_orbitals=1, i.e. 1s^2): variational theorem gives two-sided
 *    bound which can be check without any HF-specific reference number:
 *      E_exact_nonrel <= E_HF <= E_simple_product_variational
 *    The upper bound is closed-form effective-nuclear-charge result
 *    (single-parameter special case of same single-Slater-determinant ansatz
 *    family HF searches over more freely). The lower bound is known
 *    essentially-exact non-relativistic helium ground state, -2.903724 Hartree
 *    (Pekeris 1959) -> HF, being a single-determinant method, can never beat
 *     true correlated ground state.
 * 2. Beryllium (Z=4, n_orbitals=2, i.e. 1s^2 2s^2): range check against
 *    well-known non-relativistic HF-limit energy for Be, ~ -14.573 Hartree
 *    (Reference: Clementi & Roetti / standard HF tables).
 * 3. Orbital normalization: \int u_k(r)^2 dr = 1 for every converged
 *    orbital.
 * 4. Aufbau ordering: converged orbital (Fock) eigenvalues are ascending.
 */

/* NOTE: this is finite-difference/dense-diagonalization solve on a finite
 * radial grid, not a converged basis-set HF calculation. HF-limit numbers
 * depends on grid resolution (N, r_max) and should not be expected to several
 * decimal places. Tolerances below are chosen to be meaningful but
 * grid-resolution-tolerant.
 */

#include "../core/matrix.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../physics/hartree_fock.h"
#include "../physics/helium.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef RUNNING_ON_VALGRIND
#define RUNNING_ON_VALGRIND 0
#endif

#define HELIUM_EXACT_NONREL -2.903724 // Reference : Pekeris (1959), Hartree

static int check_range(double got, double lo, double hi, const char *label) {
  printf("  %s: got=%.6f  expected range=[%.6f, %.6f]\n", label, got, lo, hi);

  return (got < lo || got > hi);
}

static double integral_u2(const cvector_t *u, double dr) {
  double sum = 0.0;
  for (int i = 0; i < u->n; i++) {
    sum += u->data[i].re * u->data[i].re * dr;
  }

  return sum;
}

// Test 1: helium (Z=2, n_orbitals=1) bounded between exact and simple
// variational estimate.
static int test_helium(void) {
  int N = RUNNING_ON_VALGRIND ? 40 : 160;
  double r_min = 1e-4;
  double r_max = RUNNING_ON_VALGRIND ? 6.0 : 12.0;
  double *r = linspace(r_min, r_max, N);
  double dr = r[1] - r[0];
  int max_iter = RUNNING_ON_VALGRIND ? 50 : 200;

  hf_result_t *res =
      hartree_fock_atom_s_orbitals(r, N, 2.0, 1, 0.4, 1e-8, max_iter);

  int fail = 0;
  if (!res) {
    printf("  hartree_fock_atom_s_orbitals returned NULL\n");

    free(r);

    return 1;
  }

  printf("  converged=%d in %d iterations\n", res->converged, res->iterations);
  fail |= !res->converged;

  double e_simple_variational = helium_ground_state_energy_analytic(2.0);
  double lo = HELIUM_EXACT_NONREL - (RUNNING_ON_VALGRIND ? 0.1 : 0.01);
  double hi = e_simple_variational + (RUNNING_ON_VALGRIND ? 0.1 : 0.01);
  fail |= check_range(res->total_energy, lo, hi,
                      "He total energy (variational sandwich bound)");

  fail |= check_range(integral_u2(res->orbitals[0], dr), 0.999, 1.001,
                      "He 1s orbital normalization");

  printf("  He 1s orbital energy: %.6f Hartree\n", res->orbital_energies[0]);

  hf_result_free(res);
  free(r);

  return fail;
}

// Test 2: beryllium (Z=4, n_orbitals=2), sanity range around known HF-limit
// energy (~ -14.573 Hartree).
static int test_beryllium(void) {
  int N = RUNNING_ON_VALGRIND ? 40 : 160;
  double r_min = 1e-4;
  double r_max = RUNNING_ON_VALGRIND ? 6.0 : 10.0;
  double *r = linspace(r_min, r_max, N);
  double dr = r[1] - r[0];
  int max_iter = RUNNING_ON_VALGRIND ? 50 : 300;

  hf_result_t *res =
      hartree_fock_atom_s_orbitals(r, N, 4.0, 2, 0.4, 1e-8, max_iter);

  int fail = 0;
  if (!res) {
    printf("  hartree_fock_atom_s_orbitals returned NULL\n");

    free(r);

    return 1;
  }

  printf("  converged=%d in %d iterations\n", res->converged, res->iterations);
  fail |= !res->converged;

  double lo = RUNNING_ON_VALGRIND ? -16.0 : -15.2;
  double hi = RUNNING_ON_VALGRIND ? -12.0 : -14.0;
  fail |=
      check_range(res->total_energy, lo, hi,
                  "Be total energy (sanity range around HF limit ~-14.573)");

  for (int k = 0; k < 2; k++) {
    char label[32];
    snprintf(label, sizeof label, "Be orbital %d normalization", k);

    fail |= check_range(integral_u2(res->orbitals[k], dr), 0.999, 1.001, label);
  }

  fail |= (res->orbital_energies[0] >= res->orbital_energies[1]);
  printf("  Be orbital energies: eps_1s=%.6f  eps_2s=%.6f (Aufbau: 1s below "
         "2s)\n",
         res->orbital_energies[0], res->orbital_energies[1]);

  hf_result_free(res);
  free(r);

  return fail;
}

int main(void) {
  int failed = 0;

  printf("Helium (Z=2, 1s^2): E_exact <= E_HF <= E_simple_variational:\n");
  failed += test_helium();

  printf("Beryllium (Z=4, 1s^2 2s^2): sanity range + Aufbau ordering:\n");
  failed += test_beryllium();

  if (failed) {
    printf("FAILED (%d)\n", failed);

    return 1;
  }
  printf("PASS\n");

  return 0;
}
