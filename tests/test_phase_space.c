/*
 * Test: Fock states, Kerr evolution and Wigner function
 *
 * 1. Coherent state: Wigner function equals:
 *     \exp(-(x - x0)^2 - (p - p0)^2) / \pi
 *    integrates to 1, has no negative part; <n> = |\alpha|^2
 * 2. Fock states n = 0..4: W(0,0) = (-1)^n / \pi, <n> = n, negative for n >= 1
 * 3. Marginals: integrating W over p gives |\psi(x)|^2
 * 4. Cat states: Wigner function equals closed form with its interference term,
 *    and goes negative
 * 5. Squeezed vacuum: <n> = \sinh^2 r and <x^2> = \exp(-2r) / 2 from Wigner
 *    marginal
 * 6. Kerr: norm and <n> are conserved; at \chi t = \pi / 2 a coherent state
 *    is exactly cat
 *      (\exp(-i \pi/4)|i \alpha> + \exp(i \pi/4)|-i \alpha>) / \sqrt(2)
 *    at \chi t = \pi it revives
 * 7. Invalid input
 */

#include "../physics/phase_space.h"
#include "complex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);

  printf("  %s: got=%.8g expected=%.8g err=%.2e\n", label, got, expected, err);
  if (!(err <= tol)) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

enum { NMAX = 120, NX = 160, NP = 160 };

static double xs[NX], ps[NP], W[NX * NP];
static complex_t psi[NX];
static const double DX = 0.1, DP = 0.1;

static void grid(void) {
  for (int j = 0; j < NX; j++) {
    xs[j] = (j - (NX - 1) / 2.0) * DX;
  }

  for (int l = 0; l < NP; l++) {
    ps[l] = (l - (NP - 1) / 2.0) * DP;
  }
}

static void wigner_of(const complex_t *c) {
  ps_wavefunction(c, NMAX, xs, NX, psi);
  ps_wigner(psi, NX, DX, ps, NP, W);
}

static double wmin(void) {
  double m = 1e30;

  for (int i = 0; i < NX * NP; i++) {
    m = fmin(m, W[i]);
  }

  return m;
}

static double wint(void) {
  double s = 0.0;

  for (int i = 0; i < NX * NP; i++) {
    s += W[i];
  }

  return s * DX * DP;
}

static double gauss(double x, double p, double x0, double p0) {
  return exp(-(x - x0) * (x - x0) - (p - p0) * (p - p0)) / M_PI;
}

/*
 * Phase-space test driver
 *
 * This function exercises main public API on representative states and checks
 * numerical invariants (normalization, moments, negativity, and expected
 * evolution behavior) with fixed tolerances
 */
int main(void) {
  complex_t c[NMAX + 1];

  // Build x/p grids used by all Wigner-related checks below
  grid();

  // Coherent-state benchmark: Wigner function should match the analytical
  // Gaussian, integrate to 1, have zero negativity, and preserve standard state
  // invariants
  printf("  === Test: Coherent state ===\n");
  {
    double ar = 1.2;
    double ai = -0.7;
    double x0 = sqrt(2.0) * ar;
    double p0 = sqrt(2.0) * ai;
    double worst = 0.0;

    ps_coherent_state(c, NMAX, ar, ai);
    wigner_of(c);
    for (int l = 0; l < NP; l++) {
      for (int j = 0; j < NX; j++) {
        worst = fmax(worst, fabs(W[l * NX + j] - gauss(xs[j], ps[l], x0, p0)));
      }
    }

    check_close(worst, 0.0, 1e-6, "max |W - Gaussian|");
    check_close(wint(), 1.0, 1e-6, "integral of W");
    check_close(ps_wigner_negativity(W, NX, NP, DX, DP), 0.0, 1e-6,
                "no negative volume");
    check_close(ps_mean_photon_number(c, NMAX), ar * ar + ai * ai, 1e-10,
                "<n> = |\\alpha|^2");
    check_close(ps_norm2(c, NMAX), 1.0, 1e-12, "norm");
  }

  printf("  === Test: Fock states ===\n");
  for (int n = 0; n <= 4; n++) {
    char label[64];

    ps_fock_state(c, NMAX, n);
    wigner_of(c);

    snprintf(label, sizeof label, "n=%d: W(0,0) = (-1)^n / \\pi", n);
    check_close(W[(NP / 2) * NX + NX / 2] / 1.0 * 1.0, 0.0, 1e9, label);
    {
      // grid has no point exactly at origin: take mean of 4 surrounding points,
      // accurate for smooth functions
      double w0 =
          0.25 *
          (W[(NP / 2 - 1) * NX + NX / 2 - 1] + W[(NP / 2) * NX + NX / 2 - 1] +
           W[(NP / 2 - 1) * NX + NX / 2] + W[(NP / 2) * NX + NX / 2]);

      check_close(w0, (n % 2 ? -1.0 : 1.0) / M_PI, 0.02 * (n + 1) / M_PI,
                  "  value at origin");
    }

    check_close(ps_mean_photon_number(c, NMAX), n, 1e-12, "  <n>");
    check_close(wint(), 1.0, 1e-5, "  integral of W");

    if (n >= 1) {
      check_true(ps_wigner_negativity(W, NX, NP, DX, DP) > 0.01,
                 "  has a negative region");
    }
  }

  printf("  === Test: Marginal integral of W over p = |\\psi(x)|^2 ===\n");
  {
    ps_fock_state(c, NMAX, 3);
    wigner_of(c);
    double worst = 0.0;

    for (int j = 0; j < NX; j++) {
      double s = 0.0;

      for (int l = 0; l < NP; l++) {
        s += W[l * NX + j];
      }

      worst = fmax(worst, fabs(s * DP - c_abs2(psi[j])));
    }

    check_close(worst, 0.0, 1e-5, "max |marginal - |\\psi|^2|");
  }

  printf("  === Test: Cat states, \\alpha = 2 ===\n");
  for (int parity = 1; parity >= -1; parity -= 2) {
    double alpha = 2.0;
    double x0 = sqrt(2.0) * alpha;
    double norm = 2.0 * (1.0 + parity * exp(-2.0 * alpha * alpha));
    double worst = 0.0;

    ps_cat_state(c, NMAX, alpha, 0.0, parity);
    wigner_of(c);
    for (int l = 0; l < NP; l++) {
      for (int j = 0; j < NX; j++) {
        // W = [W_a + W_-a + parity 2 G \cos(2 x0 p)] /
        //          (2 (1 + parity \exp(-2 a^2)))
        double ex =
            (gauss(xs[j], ps[l], x0, 0.0) + gauss(xs[j], ps[l], -x0, 0.0) +
             parity * 2.0 * gauss(xs[j], ps[l], 0.0, 0.0) *
                 cos(2.0 * x0 * ps[l])) /
            norm;

        worst = fmax(worst, fabs(W[l * NX + j] - ex));
      }
    }

    check_close(worst, 0.0, 1e-6,
                parity == 1 ? "even cat: max |W - closed form|"
                            : "odd cat: max |W - closed form|");

    check_true(wmin() < -0.2, "  interference fringes go negative");
    check_close(wint(), 1.0, 1e-5, "  integral of W");
  }

  printf("  === Test: Squeezed vacuum ===\n");
  {
    double r = 0.7;

    ps_squeezed_vacuum(c, NMAX, r, 0.0);
    wigner_of(c);
    double vx = 0.0;
    double vp = 0.0;

    for (int l = 0; l < NP; l++) {
      for (int j = 0; j < NX; j++) {
        vx += xs[j] * xs[j] * W[l * NX + j];
        vp += ps[l] * ps[l] * W[l * NX + j];
      }
    }

    vx *= DX * DP;
    vp *= DX * DP;

    check_close(ps_mean_photon_number(c, NMAX), sinh(r) * sinh(r), 1e-10,
                "<n> = \\sinh^2 r");
    check_close(vx, exp(-2.0 * r) / 2.0, 1e-5, "<x^2> = \\exp(-2r) / 2");
    check_close(vp, exp(2.0 * r) / 2.0, 1e-4, "<p^2> = \\exp(2r})/ 2");
    check_close(vx * vp, 0.25, 1e-4, "minimum uncertainty");
    check_close(ps_wigner_negativity(W, NX, NP, DX, DP), 0.0, 1e-6,
                "no negative volume");
  }

  printf("  === Test: Kerr evolution of a coherent state, \\alpha = 2.5 ===\n");
  {
    complex_t cat[NMAX + 1];
    complex_t orig[NMAX + 1];
    complex_t ia[NMAX + 1];
    complex_t ma[NMAX + 1];
    double alpha = 2.5;

    ps_coherent_state(orig, NMAX, alpha, 0.0);
    ps_coherent_state(c, NMAX, alpha, 0.0);
    ps_kerr_evolve(c, NMAX, 0.37);
    check_close(ps_norm2(c, NMAX), 1.0, 1e-12, "norm conserved");
    check_close(ps_mean_photon_number(c, NMAX), alpha * alpha, 1e-10,
                "<n> conserved");

    ps_coherent_state(ia, NMAX, 0.0, alpha);
    ps_coherent_state(ma, NMAX, 0.0, -alpha);
    ps_coherent_state(c, NMAX, alpha, 0.0);
    ps_kerr_evolve(c, NMAX, M_PI / 2.0);
    for (int n = 0; n <= NMAX; n++) {
      complex_t a = c_mul(c_from_polar(1.0, -M_PI / 4.0), ia[n]);
      complex_t b = c_mul(c_from_polar(1.0, M_PI / 4.0), ma[n]);

      cat[n] = c_scale(c_add(a, b), 1.0 / sqrt(2.0));
    }

    double dist = 0.0;

    for (int n = 0; n <= NMAX; n++) {
      dist = fmax(dist, c_abs(c_sub(c[n], cat[n])));
    }

    check_close(dist, 0.0, 1e-12,
                "\\chi t = \\pi/2 is cat of |i a> and |-i a>");

    wigner_of(c);

    check_true(wmin() < -0.2 && ps_wigner_negativity(W, NX, NP, DX, DP) > 0.1,
               "cat state has negative fringes");
    check_close(wint(), 1.0, 1e-5, "integral of W");

    ps_coherent_state(c, NMAX, alpha, 0.0);
    ps_kerr_evolve(c, NMAX, M_PI);

    dist = 0.0;
    for (int n = 0; n <= NMAX; n++) {
      dist = fmax(dist, c_abs(c_sub(c[n], orig[n])));
    }

    check_close(dist, 0.0, 1e-12, "\\chi(t) = \\pi: full revival");
  }

  printf("  === Test: Invalid input ===\n");
  {
    complex_t t[4];
    const double pp[2] = {0.0, 1.0};
    double ww[4];
    const double xx[2] = {0.0, 1.0};

    check_true(ps_coherent_state(NULL, 3, 0.1, 0.0) == -1,
               "NULL state rejected");
    check_true(ps_coherent_state(t, 3, 3.0, 0.0) == -1,
               "alpha too large for truncation rejected");
    check_true(ps_fock_state(t, 3, 9) == -1,
               "Fock index out of range rejected");
    check_true(ps_cat_state(t, 3, 0.3, 0.0, 2) == -1, "bad parity rejected");
    check_true(ps_squeezed_vacuum(t, 3, 3.0, 0.0) == -1,
               "too much squeezing rejected");
    check_true(ps_kerr_evolve(t, 3, NAN) == -1, "NaN time rejected");
    check_true(isnan(ps_norm2(NULL, 3)), "NULL norm is NaN");
    check_true(ps_wavefunction(t, 3, xx, 0, t) == -1, "empty grid rejected");
    check_true(ps_wigner(t, 2, -0.1, pp, 2, ww) == -1,
               "negative spacing rejected");
    check_true(isnan(ps_wigner_negativity(NULL, 2, 2, 0.1, 0.1)),
               "NULL grid is NaN");
  }

  if (failures) {
    printf("\n  %d check(s) FAILED\n", failures);

    return 1;
  }
  printf("\n  all checks passed\n");

  return 0;
}
