/*
 * Test: RCSJ Josephson junction
 *
 * 1. Below critical current phase settles at asin(i) with zero voltage
 * 2. Overdamped limit: V = Q \sqrt(i^2 - 1) above critical current
 * 3. RK4 convergence (halving dt) and a conserved-energy check without damping
 *    or drive
 * 4. Underdamped hysteresis: up sweep switches near i = 1, down sweep retraps
 *    near 4 / (\pi Q)
 * 5. ac drive: Shapiro step, a plateau of constant mean voltage V = Omega over
 *    a range of dc current
 * 6. Invalid input
 */

#include "../physics/josephson.h"
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
  printf("  === Test: Supercurrent branch ===\n");
  {
    jj_params_t p = {1.0, 0.6, 0.0, 1.0};
    jj_state_t s = {0.3, 0.0, 0.0};
    double v = jj_mean_voltage(&s, &p, 200.0, 50.0, 0.02);

    check_close(s.phi, asin(0.6), 1e-6, "phase at asin(i)");
    check_close(v, 0.0, 1e-6, "zero voltage");
  }

  printf("  === Test: Overdamped limit ===\n");
  {
    for (int c = 0; c < 3; c++) {
      double i = (double[]){1.2, 1.6, 2.5}[c];
      jj_params_t p = {0.02, i, 0.0, 1.0};
      jj_state_t s = {0.0, 0.0, 0.0};

      // one rotation of phase
      double period = 2.0 * M_PI / (p.q * sqrt(i * i - 1.0));
      double v = jj_mean_voltage(&s, &p, 2.0 * period, 10.0 * period, 0.02);

      check_close(v / p.q, sqrt(i * i - 1.0), 0.01 * sqrt(i * i - 1.0) + 1e-3,
                  "V / Q = \\sqrt(i^2 - 1)");
    }
  }

  printf("  === Test: Integrator ===\n");
  {
    jj_params_t p = {2.0, 1.4, 0.5, 0.7};
    jj_state_t a = {0.2, 0.1, 0.0}, b = a, c = a;

    jj_advance(&a, &p, 20.0, 0.04);
    jj_advance(&b, &p, 20.0, 0.02);
    jj_advance(&c, &p, 20.0, 0.01);
    double e1 = fabs(a.phi - c.phi);
    double e2 = fabs(b.phi - c.phi);

    printf("  phi errors: dt=0.04 -> %.2e, dt=0.02 -> %.2e\n", e1, e2);
    check_true(e2 < e1 && e1 / (e2 + 1e-300) > 8.0,
               "fourth-order convergence (error falls by > 8 per halving)");

    // undamped, undriven pendulum: E = v^2/2 - \cos(\phi) is conserved
    jj_params_t f = {1e12, 0.0, 0.0, 0.0};
    jj_state_t s = {1.0, 0.0, 0.0};
    double e0 = -cos(s.phi);

    jj_advance(&s, &f, 100.0, 0.01);
    check_close(0.5 * s.v * s.v - cos(s.phi), e0, 1e-8,
                "energy conserved without damping");
  }

  printf("  === Test: Hysteresis ===\n");
  {
    const double q = 3.0;
    jj_params_t p = {q, 0.0, 0.0, 1.0};
    double up[41], down[41];

    jj_sweep(&p, 0.0, 2.0, 41, up); // steps of 0.05
    jj_sweep(&p, 2.0, 0.0, 41, down);
    int up_switch = -1;
    int down_stop = -1;

    for (int k = 0; k < 41; k++) {
      if (up_switch < 0 && up[k] > 0.1) {
        up_switch = k;
      }
      if (down_stop < 0 && down[k] < 0.05) {
        down_stop = k;
      }
    }

    double i_up = 0.05 * up_switch;
    double i_down = 2.0 - 0.05 * down_stop;

    printf("  switches up at i = %.2f, retraps at i = %.2f (4 / pi Q = %.3f)\n",
           i_up, i_down, 4.0 / (M_PI * q));
    check_true(fabs(i_up - 1.0) <= 0.1, "switching current near 1");
    check_true(i_down > 0.25 && i_down < 0.65 &&
                   fabs(i_down - 4.0 / (M_PI * q)) < 0.15,
               "retrapping current near 4 / (pi Q)");
    check_true(i_down < i_up - 0.3,
               "hysteresis: retrapping well below switching");

    check_close(up[0], 0.0, 1e-6, "no voltage at zero current (up)");
    check_true(
        up[40] > 1.5 && down[0] > 1.5 && fabs(down[40]) < 1e-3,
        "voltage state at i = 2 on both sweeps; down sweep ends at V = 0");
  }

  printf("  === Test: Shapiro step ===\n");
  {
    jj_params_t p = {0.7, 0.0, 1.2, 0.8};
    double v[41];

    jj_sweep(&p, 0.2, 2.2, 41, v);
    int run = 0;
    int best = 0;

    for (int k = 0; k < 41; k++) {
      if (fabs(v[k] - 0.8) < 2e-3) {
        run++;
        if (run > best) {
          best = run;
        }
      } else {
        run = 0;
      }
    }

    printf(
        "  longest run of points with V = \\Omega: %d (current spacing 0.05)\n",
        best);

    check_true(best >= 4, "plateau at V = \\Omega over a range of dc current");
    jj_params_t no = {0.7, 0.0, 0.0, 0.8};
    double w[41];

    jj_sweep(&no, 0.2, 2.2, 41, w);
    int has = 0;

    for (int k = 0; k < 41; k++) {
      has += fabs(w[k] - 0.8) < 2e-3;
    }

    check_true(has <= 1, "no plateau without ac drive");
  }

  printf("  === Test: Invalid input ===\n");
  {
    jj_params_t bad = {0.0, 0.5, 0.0, 1.0};
    jj_params_t ok = {1.0, 0.5, 0.0, 1.0};
    jj_state_t s = {0.0, 0.0, 0.0};
    double v[2];

    check_true(jj_step(&s, &bad, 0.01) == -1 && jj_step(&s, &ok, 0.0) == -1,
               "bad q / dt");
    check_true(isnan(jj_mean_voltage(&s, &ok, -1.0, 10.0, 0.01)) &&
                   isnan(jj_mean_voltage(&s, &ok, 1.0, 0.0, 0.01)),
               "bad averaging times");
    check_true(jj_sweep(&ok, 0.0, 1.0, 1, v) == -1 &&
                   jj_sweep(&bad, 0.0, 1.0, 2, v) == -1,
               "bad sweep");
  }

  printf(failures ? "\n  %d FAILED\n" : "\n  all checks passed\n", failures);

  return failures != 0;
}
