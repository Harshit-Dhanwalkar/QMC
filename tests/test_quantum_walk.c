/*
 * Test: discrete-time quantum walk on a line (quantum_walk.h)
 *
 * 1. One Hadamard step from |up>: P(-1) = P(+1) = 1/2
 * 2. Unitarity: norm stays 1 over hundreds of steps and every coin angle
 * 3. Ballistic spreading: <x^2>/t^2 -> 1 - \sin(\theta) (three angles)
 * 4. Symmetry and structure: symmetric start gives a symmetric distribution,
 *    peaks near +/- t / \sqrt(2), light cone |x| <= t, low centre
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

enum { N = 2048 };

static double g_p[N];
static double g_c[N];

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

// One Hadamard step from |up>: mass splits equally to x = +/-1
static void test_one_hadamard_step(void) {
  printf("  === Test: One Hadamard step from |up> ===\n");

  // Allocate state vectors once and reuse them across all subtests
  cvector_t *up = cvector_alloc(N);
  cvector_t *down = cvector_alloc(N);

  qwalk_init(up, down, QWALK_COIN_UP);
  qwalk_step(up, down, M_PI / 4, 1);
  qwalk_probability(up, down, g_p);
  check_close(g_p[N / 2 - 1], 0.5, 1e-14, "P(-1)");
  check_close(g_p[N / 2 + 1], 0.5, 1e-14, "P(+1)");
  check_close(g_p[N / 2], 0.0, 1e-14, "P(0)");

  // Clean up
  cvector_free(up);
  cvector_free(down);
}

// Unitarity: Total probability (state norm) must remain 1
static void test_unitarity(void) {
  printf("  === Test: Unitarity (norm conservation over 600 steps) ===\n");

  const double angles[] = {0.3, M_PI / 4, 1.2};
  cvector_t *up = cvector_alloc(N);
  cvector_t *down = cvector_alloc(N);

  for (int k = 0; k < 3; k++) {
    qwalk_init(up, down, QWALK_COIN_SYM);
    qwalk_step(up, down, angles[k], 600);
    check_close(qwalk_norm(up, down), 1.0, 1e-12, "norm after 600 steps");
  }

  // Clean up
  cvector_free(up);
  cvector_free(down);
}

static void test_ballistic_spreading(void) {
  printf("  === Test: Ballistic spreading <x^2>/t^2 -> 1 - sin(theta) ===\n");

  const double angles[] = {0.3, M_PI / 4, 1.2};
  cvector_t *up = cvector_alloc(N);
  cvector_t *down = cvector_alloc(N);

  for (int k = 0; k < 3; k++) {
    qwalk_init(up, down, QWALK_COIN_SYM);
    qwalk_step(up, down, angles[k], 800);
    check_close(qwalk_variance(up, down) / (800.0 * 800.0),
                qwalk_asymptotic_variance(angles[k]), 2e-4, "variance ratio");
  }

  // Clean up
  cvector_free(up);
  cvector_free(down);
}

static void test_symmetric_structure(void) {
  printf("  === Test: Symmetric start, Hadamard, t = 400 ===\n");

  cvector_t *up = cvector_alloc(N);
  cvector_t *down = cvector_alloc(N);

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
  check_close(outside, 0.0, 1e-14, "nothing outside light cone");
  check_close(peak - N / 2, 400.0 / sqrt(2.0), 8.0, "front near t/\\sqrt(2)");
  check_true(g_p[N / 2] < 0.2 * g_p[peak], "front is much higher than centre");

  // Clean up
  cvector_free(up);
  cvector_free(down);
}

static void test_biased_start(void) {
  printf("  === Test: biased start |up> ===\n");

  cvector_t *up = cvector_alloc(N);
  cvector_t *down = cvector_alloc(N);

  qwalk_init(up, down, QWALK_COIN_UP);
  qwalk_step(up, down, M_PI / 4, 400);
  check_close(qwalk_mean(up, down) / 400.0, 1.0 - 1.0 / sqrt(2.0), 2e-3,
              "drift per step");

  // Clean up
  cvector_free(up);
  cvector_free(down);
}

static void test_classical_walk(void) {
  printf("  === Test: classical walk ===\n");

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

  cvector_t *up = cvector_alloc(N);
  cvector_t *down = cvector_alloc(N);

  qwalk_init(up, down, QWALK_COIN_SYM);
  qwalk_step(up, down, M_PI / 4, 400);
  check_true(sqrt(qwalk_variance(up, down)) > 5.0 * sqrt(var),
             "quantum width > 5x classical width");

  // Clean up
  cvector_free(up);
  cvector_free(down);
}

static void test_invalid_input(void) {
  printf("  === Test: invalid input ===\n");

  cvector_t *up = cvector_alloc(N);
  cvector_t *down = cvector_alloc(N);

  check_true(qwalk_init(NULL, down, 0) == -1, "NULL");
  check_true(qwalk_init(up, down, 3) == -1, "bad coin");
  check_true(qwalk_step(up, down, NAN, 1) == -1, "NaN angle");
  check_true(qwalk_step(up, down, 0.5, 0) == -1, "steps = 0");
  check_true(qwalk_classical(g_c, N, N / 2) == -1, "t too large");
  check_true(isnan(qwalk_norm(NULL, NULL)), "norm NULL");

  cvector_free(up);
  cvector_free(down);
}

// Validation for discrete-time quantum walk behavior
// This single test entrypoint runs a sequence of physics/statistics checks
// against walk implementation and reports aggregated pass/fail status through
// global `failures` counter
// NOTE:
// - `up` and `down` hold two coin components over 1D lattice
// - `g_p` is reused as a probability workspace; `g_c` is reused for
//   classical-reference probabilities
int main(void) {
  printf(" > Quantum walk on a line tests:\n");

  test_one_hadamard_step();
  test_unitarity();
  test_ballistic_spreading();
  test_symmetric_structure();
  test_biased_start();
  test_classical_walk();
  test_invalid_input();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll quantum walk checks passed\n");

  return 0;
}
