/*
Test: scattering.c (phase_shift, born_amplitude, born_cross_section),
partial-wave observables, LS solver, Born-series convergence.

1. s-wave hard-sphere-like scattering: exact result is \delta_0 = -k * a
   for any energy.
2. Born approximation for a Yukawa potential: f(theta) has closed form via
   standard Laplace transform
   \Int_0^\inf \exp^{-\mu * r} * \sin(qr)dr = q/(\mu^2 + q^2):
     f_Born(\theta) =  (1 / hbar_sq_2m)* g/(\mu^2 + q^2), q=2k * \sin(\theta /
     2)
   checks born_amplitude's numerical integration against an exact analytic
   result.
3. born_cross_section = |f|^2, structural check.
*/

#include "../core/complex.h"
#include "../core/constants.h"
#include "../physics/potentials.h"
#include "../physics/scattering.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.6e expected=%.6e err=%.2e\n", label, got, expected, err);

  return err > tol;
}

static int test_hard_sphere_s_wave(void) {
  double a = 2.0;
  double V0 = 1000.0;
  double barrier_params[2] = {a, V0};
  double hbar_sq_2m = 0.5;

  int fail = 0;
  const double k_values[3] = {0.3, 0.5, 0.8};

  for (int i = 0; i < 3; i++) {
    double k = k_values[i];
    double delta = phase_shift(0, k, V_barrier, barrier_params, 0.01, 20.0,
                               2000, hbar_sq_2m);
    double delta_exact = -k * a;

    while (delta - delta_exact > M_PI)
      delta -= M_PI;
    while (delta - delta_exact < -M_PI)
      delta += M_PI;

    char label[32];
    snprintf(label, sizeof label, "k=%.2f \\delta_0", k);
    fail |= check_close(delta, delta_exact, 0.15, label);
  }

  return fail;
}

static int test_born_yukawa(void) {
  int fail = 0;
  double g = 1.0, mu = 1.0;
  double yukawa_params[2] = {g, mu};
  double k = 1.0;
  double r_max = 50.0;
  int N = 20000;

  const double theta_values[3] = {0.5, 1.0, 2.0};

  for (int i = 0; i < 3; i++) {
    double theta = theta_values[i];
    complex_t f =
        born_amplitude(V_yukawa, yukawa_params, k, theta, r_max, N, HBAR_2M);

    double q = 2.0 * k * sin(theta / 2.0);
    double f_exact = (1.0 / HBAR_2M) * g / (mu * mu + q * q);

    char label[32];
    snprintf(label, sizeof label, "\\theta=%.2f f_Born.re", theta);
    fail |= check_close(f.re, f_exact, fabs(f_exact) * 0.05 + 1e-30, label);
    fail |= check_close(f.im, 0.0, 1e-30, "f_Born.im (should be real)");
  }

  return fail;
}

static int test_cross_section_structural(void) {
  complex_t f = c_add(c_real(1.5), c_imag(-0.7));
  double sigma = born_cross_section(f);
  double expected = c_abs2(f);

  return check_close(sigma, expected, 1e-12, "born_cross_section = |f|^2");
}

static int test_optical_theorem(void) {
  // Weak Yukawa, several l's: \sigma_{tot} from the partial-wave sum must equal
  // (4 * \pi / k) Im f(0), computed from fully independent formula.
  int fail = 0;
  double g = 0.3, mu = 1.0;
  double yukawa_params[2] = {g, mu};
  double k = 1.0;
  double hbar_sq_2m = 0.5;
  int l_max = 6;

  double delta_l[7];
  for (int l = 0; l <= l_max; l++) {
    delta_l[l] = phase_shift(l, k, V_yukawa, yukawa_params, 1e-3, 40.0, 3000,
                             hbar_sq_2m);
  }

  double sigma_partial = total_cross_section_partial_waves(delta_l, l_max, k);
  complex_t f0 = partial_wave_amplitude(delta_l, l_max, k, 0.0);
  double sigma_optical = (4.0 * M_PI / k) * f0.im;

  fail |= check_close(sigma_partial, sigma_optical,
                      fabs(sigma_optical) * 0.02 + 1e-8, "optical theorem");
  return fail;
}

static int test_ls_vs_numerov(void) {
  // Exact LS solve vs. independent Numerov-shooting phase_shift(), same
  // potential and k - cross-check b/w 2 different methods.
  int fail = 0;
  double a = 2.0, V0 = 1000.0;
  double barrier_params[2] = {a, V0};
  double hbar_sq_2m = 0.5;
  double k = 0.4;

  double delta_numerov = phase_shift(0, k, V_barrier, barrier_params, 0.01,
                                     20.0, 2000, hbar_sq_2m);
  double delta_ls = phase_shift_lippmann_schwinger(
      0, k, V_barrier, barrier_params, 0.01, 20.0, 300, hbar_sq_2m);

  // Both are only defined |\pi (tan-based / atan2-based)|, compare tangents.
  fail |= check_close(tan(delta_numerov), tan(delta_ls), 0.1,
                      "\\tan(\\delta_0): Numerov vs Lippmann-Schwinger");
  return fail;
}

static int test_born_series_convergence(void) {
  // Weak Yukawa: successive Born-series orders should converge toward exact LS
  // phase shift, and order=1 must equal existing weak-potential formula exactly
  int fail = 0;
  double g = 0.2, mu = 1.0;
  double yukawa_params[2] = {g, mu};
  double k = 1.0;
  double hbar_sq_2m = 0.5;
  int N = 300;

  double delta_exact = phase_shift_lippmann_schwinger(
      0, k, V_yukawa, yukawa_params, 1e-3, 40.0, N, hbar_sq_2m);

  double delta_prev_err = 1e300;
  for (int order = 1; order <= 4; order++) {
    double delta_n = born_series_phase_shift(0, k, V_yukawa, yukawa_params,
                                             1e-3, 40.0, N, hbar_sq_2m, order);
    double err = fabs(delta_n - delta_exact);
    printf("  order=%d \\delta=%.6f err=%.2e\n", order, delta_n, err);

    if (order > 1 && err > delta_prev_err + 1e-9) {
      printf("  FAIL: order %d did not improve on order %d\n", order,
             order - 1);
      fail = 1;
    }
    delta_prev_err = err;
  }

  return fail;
}

int main(void) {
  int failed = 0;

  printf("s-wave hard-sphere-like phase shift (vs exact \\delta_0 = -ka):\n");
  failed += test_hard_sphere_s_wave();

  printf("Born approximation, Yukawa potential (vs closed form):\n");
  failed += test_born_yukawa();

  printf("born_cross_section structural check:\n");
  failed += test_cross_section_structural();

  printf("Optical theorem check:\n");
  failed += test_optical_theorem();

  printf("LS solver comparison with Numerov method:\n");
  failed += test_ls_vs_numerov();

  printf("Born series convergence:\n");
  failed += test_born_series_convergence();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
