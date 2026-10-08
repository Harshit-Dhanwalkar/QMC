/*
 * Test: Landau-Zener sweeps (landau_zener.h)
 *
 * 1. Closed form P_LZ = exp(-pi Omega^2 / (2 rate)) against RK4 integration for
 *    nine (\Omega, rate) pairs, starting in lower adiabatic state
 * 2. Limits: adiabatic (slow) sweep -> 0, sudden (fast) sweep -> 1, no
 *    coupling -> exactly 1, |v| conserved
 * 3. Direction: reversing sweep gives same probability
 * 4. Stueckelberg interference: two passes oscillate between ~0 and 4 P (1 - P)
 *    as sweep amplitude (hence phase) is changed 
 * 5. Step-by-step lz_advance agrees with one lz_sweep
 * 6. Invalid input
 */

#include "../physics/landau_zener.h"
#include <math.h>
#include <stdio.h>

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

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

// Upper-level population after `passes` passes from lower state
static double run(double omega, double rate, double amp, int passes) {
  double v[3];
  double end = (passes % 2 == 1) ? amp : -amp;

  lz_lower_state(v, omega, -amp);
  if (lz_sweep(v, omega, rate, amp, 0.01, passes) != 0) {
    return NAN;
  }

  return lz_upper_population(v, omega, end);
}

int main(void) {
  printf(" > Landau-Zener sweeps\n\n");

  printf("  === Closed form vs integration ===\n");
  const double omegas[] = {0.5, 1.0, 1.5};
  const double rates[] = {0.2, 1.0, 3.0};

  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 3; b++) {
      double amp = 40.0 * fmax(1.0, omegas[a]);
      double v[3];

      lz_lower_state(v, omegas[a], -amp);
      lz_sweep(v, omegas[a], rates[b], amp, 0.01, 1);
      char label[64];

      snprintf(label, sizeof label, "Omega=%.1f rate=%.1f", omegas[a],
               rates[b]);
      check_close(lz_upper_population(v, omegas[a], amp),
                  lz_probability(omegas[a], rates[b]), 2e-4, label);
      check_close(norm3(v), 1.0, 1e-6, "  |v|");
    }
  }

  printf("  === Limits ===\n");
  check_true(run(1.0, 0.02, 30.0, 1) < 1e-6,
             "adiabatic: stays on lower level");
  check_true(run(0.1, 20.0, 30.0, 1) > 0.999,
             "sudden: jumps to upper level");
  check_close(run(0.0, 1.0, 20.0, 1), 1.0, 1e-9, "Omega = 0: exactly 1");

  printf("  === Sweep direction ===\n");
  double up = run(1.0, 1.0, 40.0, 1);
  /* lz_sweep always ramps up first; reverse experiment is run with
   * lz_advance and a negative rate */
  double v2[3];

  lz_lower_state(v2, 1.0, -40.0);
  lz_advance(v2, 1.0, -40.0, 1.0, 80.0, 0.01);
  double p_forward = lz_upper_population(v2, 1.0, 40.0);

  lz_lower_state(v2, 1.0, 40.0);
  lz_advance(v2, 1.0, 40.0, -1.0, 80.0, 0.01);
  double p_backward = lz_upper_population(v2, 1.0, -40.0);

  check_close(p_forward, up, 1e-9, "forward matches lz_sweep");
  check_close(p_backward, p_forward, 2e-4, "backward = forward");

  printf("  === Stueckelberg interference (two passes) ===\n");
  double pl = lz_probability(1.0, 1.0);
  double hi = 0.0;
  double lo = 1.0;

  for (int k = 0; k < 60; k++) {
    double p = run(1.0, 1.0, 20.0 + 0.05 * k, 2);

    hi = fmax(hi, p);
    lo = fmin(lo, p);
  }
  printf("  P_LZ = %.5f, envelope 4P(1-P) = %.5f, scan range [%.5f, %.5f]\n",
         pl, 4.0 * pl * (1.0 - pl), lo, hi);
  check_true(hi <= 4.0 * pl * (1.0 - pl) + 1e-3, "never above envelope");
  check_close(hi, 4.0 * pl * (1.0 - pl), 5e-3, "reaches envelope");
  check_true(lo < 2e-3, "reaches zero (destructive)");

  printf("  === Stepwise advance equals one sweep ===\n");
  double a1[3];
  double a2[3];

  lz_lower_state(a1, 1.0, -15.0);
  lz_lower_state(a2, 1.0, -15.0);
  lz_sweep(a1, 1.0, 0.8, 15.0, 0.01, 1);
  double dur = 2.0 * 15.0 / 0.8;

  for (int i = 0; i < 10; i++) {
    lz_advance(a2, 1.0, -15.0 + 0.8 * dur * i / 10.0, 0.8, dur / 10.0, 0.01);
  }

  check_close(a1[0], a2[0], 1e-6, "x");
  check_close(a1[2], a2[2], 1e-6, "z");

  printf("  === Invalid input ===\n");
  double bad[3] = {2.0, 0.0, 0.0};

  check_true(lz_sweep(NULL, 1.0, 1.0, 5.0, 0.01, 1) == -1, "NULL v");
  check_true(lz_sweep(bad, 1.0, 1.0, 5.0, 0.01, 1) == -1, "|v| > 1");
  check_true(lz_sweep(a1, -1.0, 1.0, 5.0, 0.01, 1) == -1, "negative Omega");
  check_true(lz_sweep(a1, 1.0, 0.0, 5.0, 0.01, 1) == -1, "rate = 0");
  check_true(lz_sweep(a1, 1.0, 1.0, 0.0, 0.01, 1) == -1, "amp = 0");
  check_true(lz_sweep(a1, 1.0, 1.0, 5.0, 0.01, 0) == -1, "passes = 0");
  check_true(isnan(lz_probability(1.0, 0.0)), "P_LZ rate = 0");
  check_true(lz_lower_state(a1, 0.0, 0.0) == -1, "b = 0");

  printf("\n  %s\n", failures ? "FAILED" : "all checks passed");

  return failures ? 1 : 0;
}
