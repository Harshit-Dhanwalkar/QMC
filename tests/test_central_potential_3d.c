/*
 * Test + demonstration: central_potential_radial_solve() validated directly
 * against two closed-form/independently-computed 3D quantum systems.
 *
 * 1. 3D isotropic harmonic oscillator: exact energies
 *     E_(n,l) = \hbar * \omega*(2n + l + 3/2), n=0,1,2,..., l=0,1,2,... - about
 *   as clean closed form as exists for a 3D central potential, and exercises
 *   solver at several l simultaneously (different centrifugal barriers).
 * 2. 3D finite spherical well (l=0): bound-state energies satisfy
 *   transcendental condition
 *    k * \cot(k * a) = -\kappa (k=\sqrt(2m * (E+V0))/\hbar
 *    \kappa = \sqrt(-2m * E) / \hbar, E<0), found here via independent
 *   bisection root-finding (NOT using central_potential_radial_solve itself)
 *   before comparing to the solver's output.
 */

#include "../core/matrix.h"
#include "../core/utils.h"
#include "../physics/central_potential.h"
#include "../physics/potentials.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.6f expected=%.6f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

// Test 1: 3D isotropic harmonic oscillator, l=0,1,2, lowest 3 states each.
static int test_3d_harmonic_oscillator(void) {
  int N = 200;
  double r_min = 0.01, r_max = 8.0;
  double *r = linspace(r_min, r_max, N);
  double hbar = 1.0, mass = 1.0, omega = 1.0;

  int fail = 0;
  for (int l = 0; l <= 2; l++) {
    eigen_t *eig =
        central_potential_radial_solve(r, N, l, hbar, mass, V_harmonic, &omega);
    if (!eig) {
      printf("  l=%d: solve failed\n", l);
      fail = 1;

      continue;
    }

    for (int n = 0; n < 3; n++) {
      char label[32];
      snprintf(label, sizeof label, "l=%d n=%d", l, n);
      fail |= check_close(eig->eigenvalues[n], 2.0 * n + l + 1.5, 0.05, label);
    }

    eigen_free(eig);
  }

  free(r);

  return fail;
}

// Bisection root-find of exact transcendental bound-state condition for l=0
// finite spherical well
static double well_condition(double E, double V0, double a, double m,
                             double hbar) {
  double k = sqrt(2.0 * m * (E + V0)) / hbar;
  double kappa = sqrt(-2.0 * m * E) / hbar;

  return k / tan(k * a) + kappa;
}

static int find_well_bound_states(double V0, double a, double m, double hbar,
                                  double *out, int max_states) {
  int found = 0;
  int npts = 4000;
  double prevE = -V0 + 1e-6;
  double prevVal = well_condition(prevE, V0, a, m, hbar);

  for (int i = 1; i < npts && found < max_states; i++) {
    double E = -V0 + (double)i / npts * V0;
    double val = well_condition(E, V0, a, m, hbar);

    if ((val > 0) != (prevVal > 0) && fabs(val) < 100 && fabs(prevVal) < 100) {
      double lo = prevE, hi = E;
      double lo_sign = well_condition(lo, V0, a, m, hbar) > 0;

      for (int b = 0; b < 100; b++) {
        double mid = 0.5 * (lo + hi);
        double vmid = well_condition(mid, V0, a, m, hbar);

        if ((vmid > 0) == lo_sign) {
          lo = mid;
        } else
          hi = mid;
      }

      out[found++] = 0.5 * (lo + hi);
    }

    prevE = E;
    prevVal = val;
  }

  return found;
}

// Test 2: 3D finite spherical well, l=0 bound states.
static int test_finite_spherical_well(void) {
  double a = 3.0, V0 = 5.0, mass = 1.0, hbar = 1.0;

  double exact[8];
  int n_exact = find_well_bound_states(V0, a, mass, hbar, exact, 8);
  printf("  found %d bound states via independent root-finding\n", n_exact);

  int N = 220;
  double r_min = 0.01, r_max = 18.0;
  double *r = linspace(r_min, r_max, N);
  double params[2] = {a, V0};
  eigen_t *eig = central_potential_radial_solve(r, N, 0, hbar, mass,
                                                V_finite_well, params);

  int fail = 0;
  if (!eig) {
    printf("  solve failed\n");
    free(r);

    return 1;
  }

  for (int i = 0; i < n_exact; i++) {
    char label[32];
    snprintf(label, sizeof label, "bound state %d", i);
    fail |= check_close(eig->eigenvalues[i], exact[i], 0.05, label);
  }

  eigen_free(eig);
  free(r);

  return fail;
}

int main(void) {
  int failed = 0;

  printf("3D isotropic harmonic oscillator, E(n,l)=2n+l+3/2 "
         "(\\hbar=m=\\omega=1):\n");
  failed += test_3d_harmonic_oscillator();

  printf("3D finite spherical well (l=0), vs independent transcendental "
         "root-find:\n");
  failed += test_finite_spherical_well();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
