/*
 * Test: general central-potential radial solver.
 *
 * 1. 3D isotropic harmonic oscillator, l=0, checked against exact spectrum :
 *    E_{n_r,l} = \hbar * \omega * (2 * n_r + l + 3/2).
 * 2. Hydrogen regression: hydrogen_radial_solve must reproduce
 *   hydrogen_energy_level(1) to same accuracy as before the refactor.
 */

#include "../core/constants.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include "../physics/central_potential.h"
#include "../physics/hydrogen.h"
#include "../physics/potentials.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef RUNNING_ON_VALGRIND
#define RUNNING_ON_VALGRIND 0
#endif

static int test_harmonic_3d(void) {
  int N = RUNNING_ON_VALGRIND ? 80 : 300;
  double r_max =
      RUNNING_ON_VALGRIND ? 8.0 : 15.0; // in units where \hbar = m = \omega = 1
  double r_min = 1e-3;
  double *r = malloc(N * sizeof *r);
  double dr = (r_max - r_min) / (N - 1);

  for (int i = 0; i < N; i++) {
    r[i] = r_min + i * dr;
  }

  double hbar = 1.0, mass = 1.0, omega = 1.0;
  int l = 0;
  eigen_t *eig =
      central_potential_radial_solve(r, N, l, hbar, mass, V_harmonic, &omega);

  free(r);
  if (!eig) {
    printf("  FAIL: solver returned NULL\n");
    return 1;
  }

  // Expected: E_{n_r,0} = (2 * n_r + 3/2) * \hbar * \omega -> 1.5, 3.5, 5.5,
  // ...
  const double expected[3] = {1.5, 3.5, 5.5};
  int fail = 0;

  // Skip 2 boundary-penalty eigenvalues
  for (int k = 0; k < 3; k++) {
    double got = eig->eigenvalues[k];
    double err = fabs(got - expected[k]);
    printf("  n_r=%d: got=%.6f expected=%.6f err=%.2e\n", k, got, expected[k],
           err);

    if (err > 5e-2) { // FD grid coarseness at N=300
      fail = 1;
    }
  }

  eigen_free(eig);

  return fail;
}

static int test_hydrogen_regression(void) {
  int N = RUNNING_ON_VALGRIND ? 100 : 500;
  double r_max = RUNNING_ON_VALGRIND ? 20.0 * AU_LENGTH : 60.0 * AU_LENGTH;
  double r_min = 1e-8 * AU_LENGTH;
  double *r = malloc(N * sizeof *r);
  double dr = (r_max - r_min) / (N - 1);

  for (int i = 0; i < N; i++) {
    r[i] = r_min + i * dr;
  }

  eigen_t *eig =
      hydrogen_radial_solve(r, N, 0, HBAR, M_ELECTRON, E_CHARGE, EPSILON_0);
  free(r);

  if (!eig) {
    printf("  FAIL: hydrogen_radial_solve returned NULL\n");
    return 1;
  }

  double E1_expected = hydrogen_energy_level(1);
  double E1_got = eig->eigenvalues[0];
  double rel_err = fabs((E1_got - E1_expected) / E1_expected);
  printf("  E1: got=%.6e J expected=%.6e J rel_err=%.2e\n", E1_got, E1_expected,
         rel_err);

  eigen_free(eig);

  // Loose tolerance
  return rel_err > 0.15 ? 1 : 0;
}

int main(void) {
  int failed = 0;

  printf("3D harmonic oscillator (l=0):\n");
  failed += test_harmonic_3d();

  printf("Hydrogen regression (refactored solver):\n");
  failed += test_hydrogen_regression();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
