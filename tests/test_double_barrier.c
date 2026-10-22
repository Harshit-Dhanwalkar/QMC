/*
 * Test: transfer-matrix transmission through square barriers
 *
 * 1. Single barrier: closed forms below top, at top (E = V) and above it
 *    (\sin^2 form)
 * 2. Matching amplitudes agree with transfer-matrix formula and conserve
 *    probability (T + R = 1)
 * 3. Free space: T = 1, dwell time L / 2k; layer order reversed gives same T
 *   (reciprocity)
 * 4. Symmetric double barrier: resonant peaks of height 1 near well levels
 *    (n \pi / w)^2; asymmetric barriers lower peak to
 *    [4 GL GR / (GL + GR)^2]-like values below 1
 * 5. Probability current is constant through structure and dwell time peaks at
 *    resonance
 * 6. Tunnelling suppression at low energy; T -> 1 far above barrier
 * 7. Invalid input
 */

#include "../physics/double_barrier.h"
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

// Highest T in [lo, hi] by golden-section refinement of a coarse scan
static double peak(int n, const double *v, const double *d, double lo,
                   double hi, double *at) {
  double best = -1.0;
  double be = lo;

  for (int i = 0; i <= 4000; i++) {
    double e = lo + (hi - lo) * i / 4000;
    double t = db_transmission(n, v, d, e);

    if (t > best) {
      best = t;
      be = e;
    }
  }

  double a = be - (hi - lo) / 4000;
  double b = be + (hi - lo) / 4000;

  for (int it = 0; it < 80; it++) {
    double c = a + 0.381966 * (b - a);
    double e = a + 0.618034 * (b - a);

    if (db_transmission(n, v, d, c) > db_transmission(n, v, d, e)) {
      b = e;
    } else {
      a = c;
    }
  }

  *at = 0.5 * (a + b);

  return db_transmission(n, v, d, *at);
}

// Regression tests for transfer-matrix double-barrier routines
// Blocks below mirror numbered scenarios listed in file header:
// - analytic single-barrier limits, amplitude/probability consistency
// - reciprocity, resonance behavior, current conservation, asymptotics
// - invalid-input handling
int main(void) {
  // Single square barrier: compare against closed forms for E<V, E=V, E>V
  printf("  === Test: Single barrier ===\n");
  {
    double v = 10.0;
    double a = 0.8;

    for (int i = 0; i < 3; i++) {
      double e = (double[]){2.0, 6.5, 9.9}[i];
      double kap = sqrt(v - e);
      double want = 1.0 / (1.0 + v * v * sinh(kap * a) * sinh(kap * a) /
                                     (4 * e * (v - e)));

      check_close(db_transmission(1, &v, &a, e), want, 1e-12, "below top");
    }

    check_close(db_transmission(1, &v, &a, v), 1.0 / (1.0 + v * a * a / 4.0),
                1e-12, "at top E = V");
    for (int i = 0; i < 3; i++) {
      double e = (double[]){10.5, 14.0, 30.0}[i];
      double q = sqrt(e - v);
      double want =
          1.0 / (1.0 + v * v * sin(q * a) * sin(q * a) / (4 * e * (e - v)));

      check_close(db_transmission(1, &v, &a, e), want, 1e-12, "above top");
    }
  }

  printf("  === Test: Matching amplitudes vs transfer matrix, probability "
         "conservation\n");
  {
    const double v[3] = {12.0, 0.0, 12.0};
    const double d[3] = {0.5, 2.0, 0.7};
    double worst = 0.0;
    double worst_u = 0.0;

    for (int i = 1; i < 200; i++) {
      double e = 0.1 * i;
      db_result_t r;

      db_solve(3, v, d, e, &r);
      worst = fmax(worst, fabs(r.transmission - db_transmission(3, v, d, e)));
      worst_u = fmax(worst_u, fabs(r.transmission + r.reflection - 1.0));
    }

    printf("  max |T(matching) - T(formula)| = %.2e, max |T + R - 1| = %.2e\n",
           worst, worst_u);

    check_true(worst < 1e-9, "matching agrees with formula");
    check_true(worst_u < 1e-9, "T + R = 1");
  }

  printf("  === Test: Free space and reciprocity ===\n");
  {
    double v0 = 0.0;
    double d0 = 3.0;
    const double v[3] = {5.0, 0.0, 9.0};
    const double d[3] = {0.4, 1.5, 0.3};
    const double vr[3] = {9.0, 0.0, 5.0};
    const double dr[3] = {0.3, 1.5, 0.4};

    check_close(db_transmission(1, &v0, &d0, 2.5), 1.0, 1e-12, "free T");
    check_close(db_dwell_time(1, &v0, &d0, 4.0), 3.0 / (2 * 2.0), 1e-9,
                "free dwell time L / 2k");
    check_close(db_transmission(3, v, d, 6.3), db_transmission(3, vr, dr, 6.3),
                1e-12, "same T with order reversed");
  }

  printf("  === Test: Symmetric double barrier resonances ===\n");
  {
    const double v[3] = {40.0, 0.0, 40.0};
    const double d[3] = {0.5, 3.0, 0.5};
    double w = 3.0;

    for (int n = 1; n <= 3; n++) {
      double well = pow(n * M_PI / w, 2);
      double at;
      double p = peak(3, v, d, well * 0.75, well * 1.05, &at);

      printf("  n=%d well level %.4f, resonance at %.4f\n", n, well, at);
      check_close(p, 1.0, 1e-6, "peak transmission 1");
      check_true(fabs(at - well) / well < 0.25 && at < well,
                 "resonance just below infinite-well level");
    }

    const double v2[3] = {40.0, 0.0, 40.0};
    const double d2[3] = {0.5, 3.0, 0.8};
    double at;
    double p = peak(3, v2, d2, 0.75, 1.05, &at);

    check_true(p < 0.9 && p > 0.0, "asymmetric barriers lower peak below 1");
  }

  printf("  === Test: Current conservation and dwell time at resonance ===\n");
  {
    double v[3] = {40.0, 0.0, 40.0};
    double d[3] = {0.5, 3.0, 0.5};
    double at;

    peak(3, v, d, 0.75, 1.05, &at);
    // current ~ Im(\psi* \psi'); check via |\psi|^2 slope-free test:
    // transmitted
    // |\psi|^2 equals T in right lead
    double xs[2];
    double p2[2];
    db_result_t r;

    db_solve(3, v, d, at, &r);
    xs[0] = 4.5;
    xs[1] = 8.0;
    db_profile(3, v, d, at, xs[0], xs[1], 2, p2);

    check_close(p2[0], r.transmission, 1e-9,
                "|\\psi|^2 in right lead = T (near x = 4.5)");
    check_close(p2[1], r.transmission, 1e-9,
                "|\\psi|^2 in right lead = T (x = 8)");
    double off = db_dwell_time(3, v, d, 0.5 * at);
    double on = db_dwell_time(3, v, d, at);

    printf("  dwell time on resonance %.3f, off %.3f\n", on, off);
    check_true(on > 20.0 * off, "dwell time peaks at resonance");

    // density in well exceeds incident density at resonance
    double mid[2];

    db_profile(3, v, d, at, 2.0, 2.0 + 1e-9, 2, mid);
    check_true(mid[0] > 5.0, "density in well builds up");
  }

  printf("  === Test: Tunnelling and classical limits ===\n");
  {
    double v = 10.0;
    double a = 1.5;

    check_true(db_transmission(1, &v, &a, 0.01) < 1e-4,
               "E << V strongly suppressed");
    check_close(db_transmission(1, &v, &a, 4000.0), 1.0, 1e-4,
                "E >> V transmits");
  }

  printf("  === Test: Invalid input ===\n");
  {
    double v = 5.0;
    double a = 1.0;
    double neg = -1.0;
    double buf[4];
    db_result_t r;

    check_true(isnan(db_transmission(1, &v, &a, 0.0)) &&
                   isnan(db_transmission(1, &v, &a, -1.0)),
               "non-positive energy");
    check_true(isnan(db_transmission(0, &v, &a, 1.0)) &&
                   isnan(db_transmission(20, &v, &a, 1.0)),
               "bad layer count");
    check_true(isnan(db_transmission(1, &v, &neg, 1.0)), "negative thickness");
    check_true(db_solve(1, &v, &a, NAN, &r) == -1 &&
                   db_profile(1, &v, &a, 1.0, 1.0, 0.0, 4, buf) == -1 &&
                   db_profile(1, &v, &a, 1.0, 0.0, 1.0, 1, buf) == -1,
               "bad solve / profile arguments");
    check_true(isnan(db_dwell_time(1, &v, &a, 0.0)), "bad dwell energy");
  }

  printf(failures ? "\n  %d FAILED\n" : "\n  all checks passed\n", failures);

  return failures != 0;
}
