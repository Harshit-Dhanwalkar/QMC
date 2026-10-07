/*
 * Test: discrete-time quantum walk on a line (quantum_walk.h)
 *
 * 1. One Hadamard step from |up>: P(-1) = P(+1) = 1/2
 * 2. Unitarity: norm stays 1 over hundreds of steps and every coin angle
 * 3. Ballistic spreading: <x^2>/t^2 -> 1 - \sin(\theta) (three angles)
 * 4. Symmetry and structure: symmetric start gives a symmetric distribution,
 *    peaks near +/- t/sqrt(2), light cone |x| <= t, low centre
 * 5. Biased start (|up>): drift <x> -> (1 - 1 / \sqrt(2)) t
 * 6. Classical walk: binomial probabilities, variance t, quantum much wider
 * 7. Invalid input
 */

#include "../physics/quantum_walk.h"
#include "vector.h"
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

enum { N = 2048 };

static double g_p[N];
static double g_c[N];

int main(void) {
  printf(" > Quantum walk on a line\n\n");

  cvector_t *up = cvector_alloc(N);
  cvector_t *down = cvector_alloc(N);
  if (!up || !down) {
    return 1;
  }

  printf("  === One Hadamard step ===\n");
  qwalk_init(up, down, QWALK_COIN_UP);
  qwalk_step(up, down, M_PI / 4, 1);
  qwalk_probability(up, down, g_p);
  check_close(g_p[N / 2 - 1], 0.5, 1e-14, "P(-1)");
  check_close(g_p[N / 2 + 1], 0.5, 1e-14, "P(+1)");
  check_close(g_p[N / 2], 0.0, 1e-14, "P(0)");

  printf("  === Unitarity ===\n");
  const double angles[] = {0.3, M_PI / 4, 1.2};

  for (int k = 0; k < 3; k++) {
    qwalk_init(up, down, QWALK_COIN_SYM);
    qwalk_step(up, down, angles[k], 600);
    check_close(qwalk_norm(up, down), 1.0, 1e-12, "norm after 600 steps");
  }

  printf("  === Ballistic spreading <x^2>/t^2 -> 1 - \\sin(\\theta) ===\n");
  for (int k = 0; k < 3; k++) {
    qwalk_init(up, down, QWALK_COIN_SYM);
    qwalk_step(up, down, angles[k], 800);
    check_close(qwalk_variance(up, down) / (800.0 * 800.0),
                qwalk_asymptotic_variance(angles[k]), 2e-4, "variance ratio");
  }

  printf("  === Symmetric start, Hadamard, t = 400 ===\n");
  qwalk_init(up, down, QWALK_COIN_SYM);
  qwalk_step(up, down, M_PI / 4, 400);
  qwalk_probability(up, down, g_p);

  double asym = 0.0;
  double outside = 0.0;
  int peak = N / 2;

  for (int j = 0; j < N; j++) {
    int x = j - N / 2;

    asym = fmax(asym, fabs(g_p[j] - g_p[N - j >= N ? 0 : N - j]));
    if (abs(x) > 400) {
      outside += g_p[j];
    }
    if (x > 0 && g_p[j] > g_p[peak]) {
      peak = j;
    }
  }

  check_close(qwalk_mean(up, down), 0.0, 1e-9, "mean");
  check_close(asym, 0.0, 1e-12, "P(x) = P(-x)");
  check_close(outside, 0.0, 1e-14, "nothing outside the light cone");
  check_close(peak - N / 2, 400.0 / sqrt(2.0), 8.0, "front near t/\\sqrt(2)");
  check_true(g_p[N / 2] < 0.2 * g_p[peak],
             "front is much higher than the centre");

  printf("  === Biased start |up> ===\n");
  qwalk_init(up, down, QWALK_COIN_UP);
  qwalk_step(up, down, M_PI / 4, 400);
  check_close(qwalk_mean(up, down) / 400.0, 1.0 - 1.0 / sqrt(2.0), 2e-3,
              "drift per step");

  printf("  === Classical walk ===\n");
  check_true(qwalk_classical(g_c, N, 4) == 0, "classical t = 4");
  check_close(g_c[N / 2], 6.0 / 16.0, 1e-14, "P(0) = 6/16");
  check_close(g_c[N / 2 + 2], 4.0 / 16.0, 1e-14, "P(2) = 4/16");
  check_close(g_c[N / 2 + 4], 1.0 / 16.0, 1e-14, "P(4) = 1/16");
  qwalk_classical(g_c, N, 400);
  double var = 0.0;
  double tot = 0.0;

  for (int j = 0; j < N; j++) {
    var += (j - N / 2.0) * (j - N / 2.0) * g_c[j];
    tot += g_c[j];
  }

  check_close(tot, 1.0, 1e-12, "classical norm");
  check_close(var, 400.0, 1e-8, "classical variance = t");
  qwalk_init(up, down, QWALK_COIN_SYM);
  qwalk_step(up, down, M_PI / 4, 400);
  check_true(sqrt(qwalk_variance(up, down)) > 5.0 * sqrt(var),
             "quantum width > 5x classical width");

  printf("  === Invalid input ===\n");
  check_true(qwalk_init(NULL, down, 0) == -1, "NULL");
  check_true(qwalk_init(up, down, 3) == -1, "bad coin");
  check_true(qwalk_step(up, down, NAN, 1) == -1, "NaN angle");
  check_true(qwalk_step(up, down, 0.5, 0) == -1, "steps = 0");
  check_true(qwalk_classical(g_c, N, N / 2) == -1, "t too large");
  check_true(isnan(qwalk_norm(NULL, NULL)), "norm NULL");

  cvector_free(up);
  cvector_free(down);

  printf("\n  %s\n", failures ? "FAILED" : "all checks passed");

  return failures ? 1 : 0;
}
