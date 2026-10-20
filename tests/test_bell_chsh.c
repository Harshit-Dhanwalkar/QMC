/*
 * Test: CHSH Bell test on a Werner state
 *
 * 1. Density-matrix correlation equals -V \cos(a - b)
 * 2. Standard angles give S = 2 sqrt 2 V (Tsirelson); grid search finds no
 *    larger value
 * 3. Violation threshold: S > 2 exactly when V > 1 / sqrt 2
 * 4. deterministic-strategy bound is 2, and hidden-angle local model
 *    never exceeds it
 * 5. Local model's correlation matches a direct Monte-Carlo of sign(\cos(a -
 *    \lambda))
 * 6. Sampled correlations converge to quantum values with 1 / \sqrt{n} error
 * 7. Invalid input
 */

#include "../physics/bell_chsh.h"
#include <math.h>
#include <stdio.h>

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

int main(void) {
  const double A = 0.0;
  const double A2 = M_PI / 2.0;
  const double B = M_PI / 4.0;
  const double B2 = 3.0 * M_PI / 4.0;

  // Validate analytical Werner-state correlation formula
  // E(a,b) = -V \cos(a-b) across a sweep of representative angles/visibilities
  // Also check two anchor points (perfect anticorrelation for V=1 on same axis,
  // and zero correlation for maximally mixed case V=0)
  printf("  === Test: Correlation from the density matrix ===\n");
  {
    double worst = 0.0;

    for (int i = 0; i < 40; i++) {
      double a = 0.37 * i;
      double b = 1.1 - 0.23 * i;
      double v = 0.1 + 0.02 * i;

      worst = fmax(worst, fabs(bell_correlation(v, a, b) + v * cos(a - b)));
    }

    check_close(worst, 0.0, 1e-14, "max |E + V cos(a - b)|");
    check_close(bell_correlation(1.0, 0.7, 0.7), -1.0, 1e-14,
                "singlet, same axis: perfect anticorrelation");
    check_close(bell_correlation(0.0, 0.3, 1.9), 0.0, 1e-14,
                "pure noise: no correlation");
  }

  // Verify CHSH behavior at canonical settings:
  // - magnitude equals 2*sqrt(2)*V at the standard optimal angles
  // - grid search reproduces the same maximum
  // - violation threshold occurs only above V > 1/\sqrt(2)
  printf("  === Test: Tsirelson bound ===\n");
  {
    check_close(fabs(bell_chsh(1.0, A, A2, B, B2)), 2.0 * M_SQRT2, 1e-14,
                "S at the standard angles, V = 1");
    check_close(fabs(bell_chsh(0.8, A, A2, B, B2)), 2.0 * M_SQRT2 * 0.8, 1e-14,
                "S at V = 0.8");
    check_close(bell_max_chsh(1.0, 16), 2.0 * M_SQRT2, 1e-12,
                "grid search maximum, V = 1");
    check_close(bell_max_chsh(0.5, 16), M_SQRT2, 1e-12,
                "grid search maximum, V = 0.5");

    check_true(fabs(bell_chsh(0.7, A, A2, B, B2)) < 2.0,
               "V = 0.7 < 1/sqrt 2: no violation");
    check_true(fabs(bell_chsh(0.72, A, A2, B, B2)) > 2.0,
               "V = 0.72 > 1/sqrt 2: violation");
  }

  // Cross-check local realistic baselines:
  // deterministic strategies must satisfy |S| <= 2 and hidden-angle local model
  // should not exceed this classical CHSH bound
  printf("  === Test: Local models ===\n");
  {
    double worst = 0.0;

    check_close(bell_deterministic_bound(), 2.0, 0.0,
                "deterministic strategies reach exactly 2");
    for (int i = 0; i < 30; i++) {
      for (int j = 0; j < 30; j++) {
        worst =
            fmax(worst, fabs(bell_local_chsh(0.2 * i, 0.9 * j, 0.11 * i + 0.4,
                                             0.5 * j - 0.3)));
      }
    }

    check_true(worst <= 2.0 + 1e-12, "hidden-angle model never exceeds 2");

    check_close(fabs(bell_local_chsh(A, A2, B, B2)), 2.0, 1e-14,
                "hidden-angle model at the standard angles");

    // direct Monte-Carlo of the model: hidden angle lambda, outputs sign(cos(a
    // - \lambda)), -sign(\cos(b - \lambda))
    int n = 400000;
    double a = 0.3;
    double b = 1.2;
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
      double lam = 2.0 * M_PI * (i + 0.5) / n;
      int ra = cos(a - lam) >= 0 ? 1 : -1;
      int rb = cos(b - lam) >= 0 ? -1 : 1;

      sum += ra * rb;
    }

    check_close(sum / n, bell_local_correlation(a, b), 1e-4,
                "sawtooth correlation = model average");
  }

  printf("  === Test: Sampling ===\n");
  {
    double sums[4];
    int n = 200000;
    double s = bell_sample(1.0, A, A2, B, B2, n, 12345u, sums);

    check_close(s, -2.0 * M_SQRT2, 5.0 * 4.0 / sqrt(n),
                "sampled S within 5 \\sigma of -2 \\sqrt{2}");

    for (int c = 0; c < 4; c++) {
      const double ang[4][2] = {{A, B}, {A, B2}, {A2, B}, {A2, B2}};

      check_close(sums[c] / n, bell_correlation(1.0, ang[c][0], ang[c][1]),
                  5.0 / sqrt(n), "  sampled E near exact");
    }

    double s2 = bell_sample(1.0, A, A2, B, B2, n, 12345u, sums);

    check_close(s2, s, 0.0, "same seed reproduces the sample");
  }

  printf("  === Test: Invalid input ===\n");
  {
    double sums[4];

    check_true(isnan(bell_correlation(1.2, 0.0, 0.0)), "V > 1 rejected");
    check_true(isnan(bell_correlation(-0.1, 0.0, 0.0)), "V < 0 rejected");
    check_true(isnan(bell_correlation(0.5, NAN, 0.0)), "NaN angle rejected");
    check_true(isnan(bell_max_chsh(0.5, 1)), "grid of 1 rejected");
    check_true(isnan(bell_sample(0.5, 0, 1, 2, 3, 0, 1u, sums)),
               "zero pairs rejected");
    check_true(isnan(bell_sample(0.5, 0, 1, 2, 3, 10, 1u, NULL)),
               "NULL output rejected");
  }

  if (failures) {
    printf("\n  %d check(s) FAILED\n", failures);

    return 1;
  }
  printf("\n  all checks passed\n");

  return 0;
}
