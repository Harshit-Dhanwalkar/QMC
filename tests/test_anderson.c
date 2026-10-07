/*
 * Test: 1D Anderson localisation (anderson.h)
 *
 * 1. Clean chain: Bessel-function spreading, <x^2> = 2 t^2 \tau^2 exactly,
 *    norm and energy conserved
 * 2. Disorder generator: determinism, range, mean and variance W^2/12
 * 3. Strong disorder: packet stops spreading (width saturates), IPR
 *    stays large, while a clean chain keeps spreading
 * 4. Lyapunov exponent: transfer matrices vs weak-disorder formula
 *    W^2 / (24 (4 t^2 - E^2))
 * 5. Invalid input
 */

#include "../physics/anderson.h"
#include "complex.h"
#include "vector.h"
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

enum { N = 801, MID = 400, BIG = 400000 };

static double g_eps[N];
static double g_big[BIG];

static void start_on_site(cvector_t *psi, int site) {
  cvector_fill(psi, c_new(0.0, 0.0));
  psi->data[site] = c_new(1.0, 0.0);
}

int main(void) {
  printf(" > 1D Anderson localisation\n\n");
  cvector_t *psi = cvector_alloc(N);

  if (!psi) {
    return 1;
  }

  printf("  === Clean chain: ballistic spreading ===\n");
  for (int j = 0; j < N; j++) {
    g_eps[j] = 0.0;
  }

  start_on_site(psi, MID);

  // t = 1, tau = 20 (2000 steps of 0.01); packet front at |j| ~ 2 t tau = 40
  check_true(anderson_evolve(psi, g_eps, 1.0, 0.01, 2000) == 0, "evolve");
  check_close(anderson_norm(psi), 1.0, 1e-9, "norm");
  check_close(anderson_width(psi), sqrt(2.0) * 20.0, 1e-5,
              "width = \\sqrt(2) t \\tau");
  check_close(anderson_energy(psi, g_eps, 1.0), 0.0, 1e-9, "energy");

  // hopping scales time: t = 0.5, \tau = 40 gives same t \tau
  start_on_site(psi, MID);
  anderson_evolve(psi, g_eps, 0.5, 0.01, 4000);
  check_close(anderson_width(psi), sqrt(2.0) * 0.5 * 40.0, 1e-5,
              "width with t = 0.5");

  printf("  === Disorder generator ===\n");
  double a[8];
  double b[8];

  anderson_disorder(a, 8, 3.0, 42);
  anderson_disorder(b, 8, 3.0, 42);
  int same = 1;

  for (int i = 0; i < 8; i++) {
    same &= (a[i] == b[i]);
  }

  check_true(same, "same seed, same sequence");
  anderson_disorder(b, 8, 3.0, 43);
  check_true(a[0] != b[0], "different seed differs");
  anderson_disorder(g_big, BIG, 3.0, 7);

  double lo = 1e9;
  double hi = -1e9;
  double mean = 0.0;
  double var = 0.0;

  for (int i = 0; i < BIG; i++) {
    lo = fmin(lo, g_big[i]);
    hi = fmax(hi, g_big[i]);
    mean += g_big[i];
  }

  mean /= BIG;
  for (int i = 0; i < BIG; i++) {
    var += (g_big[i] - mean) * (g_big[i] - mean);
  }

  var /= BIG;
  check_true(lo >= -1.5 && hi < 1.5, "values in [-W/2, W/2)");
  check_close(mean, 0.0, 0.01, "mean");
  check_close(var, 9.0 / 12.0, 0.01, "variance W^2/12");

  printf("  === Strong disorder localises ===\n");
  anderson_disorder(g_eps, N, 6.0, 2024);
  start_on_site(psi, MID);
  anderson_evolve(psi, g_eps, 1.0, 0.01, 2000); // \tau = 20
  double w20 = anderson_width(psi);
  double e0 = anderson_energy(psi, g_eps, 1.0);

  anderson_evolve(psi, g_eps, 1.0, 0.01, 4000); // \tau = 60
  double w60 = anderson_width(psi);

  printf("  width(20) = %.4f, width(60) = %.4f\n", w20, w60);
  check_true(w60 < 12.0, "width stays small (clean chain would be ~85)");
  check_true(fabs(w60 - w20) < 0.25 * w20 + 1.0, "width saturated");
  check_true(anderson_ipr(psi) > 0.05, "IPR large: few sites occupied");
  check_close(anderson_norm(psi), 1.0, 1e-6, "norm (RK4 drift)");
  check_close(anderson_energy(psi, g_eps, 1.0), e0, 1e-6, "energy conserved");

  printf("  === Lyapunov exponent vs weak-disorder formula ===\n");
  anderson_disorder(g_big, BIG, 1.0, 99);
  double lam = anderson_lyapunov(g_big, BIG, 1.0, 1.0);
  double ref = anderson_lyapunov_weak(1.0, 1.0, 1.0);

  printf("  W = 1, E = 1: measured %.5f, formula %.5f\n", lam, ref);
  check_true(fabs(lam / ref - 1.0) < 0.12, "within 12% (statistical + O(W^4))");
  anderson_disorder(g_big, BIG, 0.0, 1);
  check_close(anderson_lyapunov(g_big, BIG, 1.0, 1.0), 0.0, 1e-3,
              "clean chain: no growth inside band");
  check_true(isnan(anderson_lyapunov_weak(1.0, 1.0, 2.5)),
             "outside band: NaN");

  printf("  === Invalid input ===\n");
  check_true(anderson_disorder(NULL, 4, 1.0, 1) == -1, "NULL eps");
  check_true(anderson_disorder(a, 8, -1.0, 1) == -1, "negative W");
  check_true(anderson_evolve(NULL, g_eps, 1.0, 0.01, 1) == -1, "NULL psi");
  check_true(anderson_evolve(psi, g_eps, 0.0, 0.01, 1) == -1, "hop = 0");
  check_true(anderson_evolve(psi, g_eps, 1.0, 0.0, 1) == -1, "dt = 0");
  check_true(anderson_evolve(psi, g_eps, 1.0, 0.01, 0) == -1, "steps = 0");
  check_true(isnan(anderson_width(NULL)), "width NULL");
  cvector_free(psi);

  printf("\n  %s\n", failures ? "FAILED" : "all checks passed");

  return failures ? 1 : 0;
}
