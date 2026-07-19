/*
Test: scattering.c (phase_shift, born_amplitude, born_cross_section).

1. s-wave hard-sphere-like scattering: exact result is delta_0 = -k*a
   for ANY energy. No hard-wall potential exists in potentials.c, so
   this is approximated via V_barrier with a large (but finite) V0 -
   flagged as an approximation, not an exact match, with a loose
   tolerance reflecting that.
2. Born approximation for a Yukawa potential: f(theta) has a clean
   closed form via the standard Laplace transform
   \Int_0^\inf \exp^{-\mu * r} * \sin(qr)dr = q/(\mu^2 + q^2):
     f_Born(\theta) = (2 * M_ELECTRON/ \hbar^2) * g/(\mu^2 + q^2), q=2k *
     \sin(\theta/2) checks born_amplitude's numerical
   integration against an exact analytic result.
3. born_cross_section = |f|^2, a trivial structural check.

// HACK: born_amplitude hardcodes M_ELECTRON and HBAR internally rather than
taking mass/hbar as parameters

// FIX: since fixing signature requires scattering.h
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
  double k_values[3] = {0.3, 0.5, 0.8};
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
    snprintf(label, sizeof label, "k=%.2f delta_0", k);
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

  double theta_values[3] = {0.5, 1.0, 2.0};
  for (int i = 0; i < 3; i++) {
    double theta = theta_values[i];
    complex_t f = born_amplitude(V_yukawa, yukawa_params, k, theta, r_max, N);

    double q = 2.0 * k * sin(theta / 2.0);
    double f_exact = (2.0 * M_ELECTRON / (HBAR * HBAR)) * g / (mu * mu + q * q);

    char label[32];
    snprintf(label, sizeof label, "theta=%.2f f_Born.re", theta);
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

int main(void) {
  int failed = 0;

  printf("s-wave hard-sphere-like phase shift (vs exact delta_0=-ka):\n");
  failed += test_hard_sphere_s_wave();

  printf("Born approximation, Yukawa potential (vs closed form):\n");
  failed += test_born_yukawa();

  printf("born_cross_section structural check:\n");
  failed += test_cross_section_structural();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
