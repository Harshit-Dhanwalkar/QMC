/*
 * Test: Quantum Zeno effect
 *
 * 1. N = 1 pulsed measurement at T = \pi / \Omega: survival 0; N = 2: 1/4;
 *    closed form \cos^{2N}
 * 2. Survival rises monotonically to 1 as N grows and follows
 *    \exp(-(\Omega T)^2  / 4N) at large N
 * 3. Bloch equations: \gamma = 0 gives \cos^2(\Omega t / 2); long-time limit
 *    1/2
 * 4. Slowest rate: matches numerical decay of envelope; maximum \Omega at
 *    \gamma = 2 \Omega; \Omega^2 / \gamma at strong measurement
 * 5. Monte-Carlo trajectories agree with survival probability within statistics
 * 6. Invalid input
 */

#include "../physics/zeno.h"
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

// Integration-style test driver for Zeno module
// Verifies closed-form pulsed survival, time-resolved pulsed survival and
// continuous-measurement behavior
// Each failed check increments global `failures` counter and is summarized at
// process exit
int main(void) {
  const double omega = 1.0;
  const double T = M_PI;

  // Pulsed-measurement baseline checks and asymptotic Zeno behavior
  printf("  === Test: Pulsed measurements ===\n");
  {
    // Exact/reference values for small N at T = \pi/\Omega
    check_close(zeno_survival_pulsed(omega, T, 1), 0.0, 1e-30,
                "N = 1 at T = \\pi/\\Omega: complete flip");
    check_close(zeno_survival_pulsed(omega, T, 2), 0.25, 1e-14,
                "N = 2: \\cos^4(\\pi/4) = 1/4");
    check_close(zeno_survival_pulsed(omega, T, 4), pow(cos(M_PI / 8.0), 8.0),
                1e-14, "N = 4");

    int mono = 1;

    // For fixed total time T, increasing measurement count should not reduce
    // survival
    for (int n = 2; n < 400; n++) {
      if (zeno_survival_pulsed(omega, T, n + 1) <
          zeno_survival_pulsed(omega, T, n)) {
        mono = 0;
      }
    }

    check_true(mono, "survival increases with N");

    check_close(zeno_survival_pulsed(omega, T, 1000), 1.0, 5e-3,
                "N = 1000: nearly frozen");
    check_close(zeno_survival_pulsed(omega, T, 1000),
                exp(-omega * omega * T * T / 4000.0), 1e-5,
                "large-N law \\exp(-(\\Omega T)^2 / 4N)");
  }

  // Time-resolved pulsed survival should match boundary and mid-interval
  // formulas
  printf("  === Test: Survival during pulsed run ===\n");
  {
    check_close(zeno_survival_pulsed_at(omega, T, 5, 0.0), 1.0, 1e-15, "t = 0");
    check_close(zeno_survival_pulsed_at(omega, T, 5, T),
                zeno_survival_pulsed(omega, T, 5), 1e-12,
                "t = T equals total survival");
    check_close(zeno_survival_pulsed_at(omega, T, 1, 0.5 * T),
                pow(cos(omega * 0.5 * T / 2.0), 2.0), 1e-14,
                "N = 1 mid-run is free Rabi evolution");
    check_close(zeno_survival_pulsed_at(omega, T, 4, 0.3 * T),
                pow(cos(omega * T / 8.0), 2.0) *
                    pow(cos(omega * 0.05 * T / 2.0), 2.0),
                1e-12, "mid-interval after one measurement");
  }

  // Continuous-measurement model checks (Bloch-equation implementation)
  printf("  === Test: Continuous measurement: Bloch equations ===\n");
  {
    double worst = 0.0;

    for (int i = 0; i <= 20; i++) {
      double t = 0.3 * i;

      worst = fmax(worst, fabs(zeno_survival_continuous(omega, 0.0, t) -
                               pow(cos(omega * t / 2.0), 2.0)));
    }

    check_close(worst, 0.0, 1e-9, "\\gamma = 0: \\cos^2(\\Omega t / 2)");
    check_close(zeno_survival_continuous(omega, 0.7, 80.0), 0.5, 1e-6,
                "long-time limit 1/2");
  }

  printf("  === Test: Slowest decay rate ===\n");
  {
    double gammas[] = {0.4, 1.0, 2.0, 2.5, 4.0, 10.0, 40.0};

    for (unsigned i = 0; i < sizeof gammas / sizeof gammas[0]; i++) {
      double g = gammas[i];
      double r = zeno_continuous_rate(omega, g);
      // late-time envelope of (P - 1/2): ratio over a window deep in tail
      double t1 = 6.0 / fmax(r, 0.05), t2 = t1 + 1.0 / fmax(r, 0.05);
      double e1 = fabs(zeno_survival_continuous(omega, g, t1) - 0.5),
             e2 = fabs(zeno_survival_continuous(omega, g, t2) - 0.5);
      char label[80];

      if (g == 2.0 * omega) {
        snprintf(label, sizeof label,
                 "\\gamma = %.1f: continuous with overdamped branch", g);
        check_close(r, zeno_continuous_rate(omega, g + (g == 2.0 ? 1e-9 : 0.0)),
                    1e-4, label);
      } else if (g <= 2.0 * omega) {
        // underdamped: envelope exp(-gamma t / 2) is read off at times one
        // damped period apart
        double wd = sqrt(omega * omega - g * g / 4.0);
        double per = 2.0 * M_PI / wd;
        double tt1 = 3.0 * per;
        double tt2 = 4.0 * per;
        double u1 = zeno_survival_continuous(omega, g, tt1) - 0.5;
        double u2 = zeno_survival_continuous(omega, g, tt2) - 0.5;

        snprintf(label, sizeof label,
                 "\\gamma = %.1f: numerical envelope decay rate", g);
        check_close(log(fabs(u1 / u2)) / (tt2 - tt1), r, 2e-3, label);
      } else {
        snprintf(label, sizeof label,
                 "\\gamma = %.1f: numerical envelope decay rate", g);
        check_close(log(e1 / e2) / (t2 - t1), r, 2e-3, label);
      }
    }

    check_close(zeno_continuous_rate(omega, 2.0), omega, 1e-14,
                "fastest decay at \\gamma = 2 \\Omega equals \\Omega");
    check_close(zeno_continuous_rate(omega, 1000.0) * 1000.0, omega * omega,
                1e-3, "strong measurement: \\Omega^2 / \\gamma");
    check_true(zeno_continuous_rate(omega, 3.0) <
                   zeno_continuous_rate(omega, 2.0),
               "decay slows beyond \\gamma = 2 \\Omega");
  }

  printf("  === Test: Monte-Carlo trajectories ===\n");
  {
    int n = 20;
    int ntraj = 200000;
    static int alive[64];

    zeno_sample_pulsed(omega, T, n, ntraj, 99u, alive);
    check_close((double)alive[n] / ntraj, zeno_survival_pulsed(omega, T, n),
                5.0 * 0.5 / sqrt(ntraj), "N = 20 survival within 5 \\sigma");
    check_close((double)alive[1] / ntraj, pow(cos(M_PI / 40.0), 2.0),
                5.0 * 0.5 / sqrt(ntraj), "first measurement");
    check_true(alive[0] == ntraj, "all alive before first measurement");

    int mono = 1;

    for (int k = 1; k <= n; k++) {
      if (alive[k] > alive[k - 1]) {
        mono = 0;
      }
    }

    check_true(mono, "alive count never increases");
  }

  printf("  === Test: invalid input ===\n");
  {
    int a[4];

    check_true(isnan(zeno_survival_pulsed(1.0, 1.0, 0)), "N = 0 rejected");
    check_true(isnan(zeno_survival_pulsed(1.0, -1.0, 4)),
               "negative time rejected");
    check_true(isnan(zeno_survival_continuous(1.0, -0.5, 1.0)),
               "negative dephasing rejected");
    check_true(isnan(zeno_continuous_rate(NAN, 1.0)), "NaN rate rejected");
    check_true(zeno_sample_pulsed(1.0, 1.0, 3, 10, 1u, NULL) == -1,
               "NULL output rejected");
    check_true(zeno_sample_pulsed(1.0, 1.0, 0, 10, 1u, a) == -1,
               "N = 0 trajectories rejected");
  }

  if (failures) {
    printf("\n  %d check(s) FAILED\n", failures);

    return 1;
  }
  printf("\n  all checks passed\n");

  return 0;
}
