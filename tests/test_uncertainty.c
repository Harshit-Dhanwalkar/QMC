/*
 * Test: physics/uncertainty.c's compute_uncertainties and
 * compute_energy_expectation
 *
 * test_wavefn.c already covers these two functions on the quantum harmonic
 * oscillator (QHO) ground state, where the Heisenberg bound is saturated
 * (it's a minimum-uncertainty Gaussian state). This test adds a second,
 * qualitatively different analytic system - ground state of infinite square
 * well - which is smooth but non-Gaussian, so:
 *
 *   - The uncertainty product should come out clearly above \hbar/2 (not
 *     saturated)
 *   - The exact analytic value of \Delta_x, \Delta_p, and the ground-state
 *     energy are all theroretical closed forms, so this cross-checks underlying
 *     wavefunction_expect_* machinery against a second, unrelated derivation
 *
 * Also covers NULL-input safety path of compute_uncertainties
 */

#include "../physics/potentials.h"
#include "../physics/uncertainty.h"
#include "../physics/wavefn.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  if (fabs(got - expected) > tol) {
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

static void check_true(int cond, const char *msg) {
  printf("  %s: %s\n", msg, cond ? "ok" : "FAILED");
  if (!cond) {
    failures++;
  }
}

/*
 * Ground state of the infinite square well V_infinite_well (V=0 for
 * |x|<a, V=huge otherwise), natural units \hbar=m=1:
 *   \psi(x) = \sqrt(1/a) * \cos(\pi * x / (2a)),   |x| <= a
 * which vanishes at both edges (x = \pm a) and is normalized on [-a,a]:
 *   \int_{-a}^{a} (1/a) \cos^2(\pi * x / (2a)) dx = 1
 *
 * Sampling the grid over exactly [-a,a] keeps the wavefunction at (near) zero
 * at both boundaries, which is what wavefunction_expect_p/p2's FFT approach
 * needs to behave like a periodic function
 */
static wavefunction_t *make_infinite_well_ground_state(int n, double a) {
  wavefunction_t *wf = wavefunction_alloc(n);
  wf->dx = (2.0 * a) / (n - 1);

  double norm_const = sqrt(1.0 / a);
  for (int i = 0; i < n; i++) {
    double xi = -a + i * wf->dx;

    wf->x[i] = xi;
    wf->psi->data[i].re = norm_const * cos(M_PI * xi / (2.0 * a));
    wf->psi->data[i].im = 0.0;
  }

  return wf;
}

static void test_null_wavefunction_is_safe(void) {
  printf("Test: compute_uncertainties(NULL) is safe and returns all-zero\n");

  uncertainty_t u = compute_uncertainties(NULL);

  check_close(u.mean_x, 0.0, 0.0, "mean_x is zero for NULL input");
  check_close(u.mean_p, 0.0, 0.0, "mean_p is zero for NULL input");
  check_close(u.delta_x, 0.0, 0.0, "delta_x is zero for NULL input");
  check_close(u.delta_p, 0.0, 0.0, "delta_p is zero for NULL input");
  check_close(u.product, 0.0, 0.0, "product is zero for NULL input");

  // compute_energy_expectation delegates to variational_energy, which return
  // 0.0 for a NULL wavefunction/potential
  check_close(compute_energy_expectation(NULL, V_harmonic, NULL, 1.0), 0.0, 0.0,
              "compute_energy_expectation(NULL, ...) returns 0.0");
}

static void test_infinite_well_ground_state_uncertainty(void) {
  printf("Test: compute_uncertainties on the infinite-square-well ground "
         "state is not saturated (non-Gaussian state)\n");

  double a = 1.0;
  wavefunction_t *wf = make_infinite_well_ground_state(4096, a);

  // Position-space normalization sanity check (Riemann sum over dx)
  double norm = 0.0;
  for (int i = 0; i < wf->n; i++) {
    norm += wf->psi->data[i].re * wf->psi->data[i].re * wf->dx;
  }
  check_close(norm, 1.0, 1e-3, "position-space normalization holds");

  uncertainty_t u = compute_uncertainties(wf);

  /* Analytic results (natural units, hbar=m=1), well half-width a=1:
   *  <x> = 0, <p> = 0 (real, even wavefunction)
   *  Var(x) = (a^2/3) * (1 - 6/\pi^2)  =>  \delta_x = \sqrt(Var(x))
   *  E_1 = \pi^2 / (8 a^2)  =>  <p^2> = 2*E_1  =>  \delta_p = \sqrt(<p^2>)
   */
  double var_x = (a * a / 3.0) * (1.0 - 6.0 / (M_PI * M_PI));
  double delta_x_exact = sqrt(var_x);
  double E1_exact = (M_PI * M_PI) / (8.0 * a * a);
  double delta_p_exact = sqrt(2.0 * E1_exact);
  double product_exact = delta_x_exact * delta_p_exact;

  printf("  <x>=%.6f <p>=%.6f delta_x=%.6f (exact %.6f) delta_p=%.6f (exact "
         "%.6f)\n",
         u.mean_x, u.mean_p, u.delta_x, delta_x_exact, u.delta_p,
         delta_p_exact);

  check_close(u.mean_x, 0.0, 1e-6, "<x>=0 by symmetry");
  check_close(u.mean_p, 0.0, 1e-6, "<p>=0 by symmetry");
  check_close(u.delta_x, delta_x_exact, 1e-3,
              "\\delta_x matches the analytic infinite-well result");
  check_close(u.delta_p, delta_p_exact, 2e-2,
              "\\delta_p matches the analytic infinite-well result");
  check_close(u.product, product_exact, 2e-2,
              "uncertainty product matches the analytic value");

  // The physically important, qualitative check: a non-Gaussian bound state
  // must not saturate Heisenberg bound way QHO ground state  does
  check_true(u.product > 0.5 + 1e-3,
             "uncertainty product is strictly above hbar/2 (not saturated, "
             "unlike QHO ground state)");

  wavefunction_free(wf);
}

static void test_infinite_well_energy_expectation(void) {
  printf("Test: compute_energy_expectation on the infinite-well ground state "
         "matches E_1 = pi^2/(8a^2)\n");

  double a = 1.0;
  wavefunction_t *wf = make_infinite_well_ground_state(4096, a);

  double E1_exact = (M_PI * M_PI) / (8.0 * a * a);
  double E1_got = compute_energy_expectation(wf, V_infinite_well, &a, 1.0);

  printf("  E1 got=%.6f exact=%.6f\n", E1_got, E1_exact);
  check_close(E1_got, E1_exact, 2e-2,
              "energy expectation matches the exact ground-state energy");

  wavefunction_free(wf);
}

int main(void) {
  test_null_wavefunction_is_safe();
  test_infinite_well_ground_state_uncertainty();
  test_infinite_well_energy_expectation();

  if (failures > 0) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }

  printf("\nAll checks passed\n");
  return 0;
}
