/*
 * Test: Deutsch-Jozsa and Grover's search algorithms.
 *
 * 1. Deutsch-Jozsa: DJ_CONSTANT_0/DJ_CONSTANT_1 must give P(all-zero)=1.0
 *    exactly; DJ_BALANCED (single-bit and multi-bit parity oracles) must give
 *    P(all-zero)=0.0 exactly : which are exact quantum-mechanical predictions
 *    (no measurement noise involved, deterministic outcome probabilities).
 * 2. Grover's search: target's measurement probability after optimal number of
 *    iterations must match the closed-form \sin^2((2k + 1) * \theta),
 *    \sin(\theta)=1 / \sqrt(N), across several (n_qubits, target) combinations.
 *    Also checks that full probability distribution sums to 1 (normalization)
 *    and that iterating well past optimum ("over-rotation") measurably reduces
 *    the success probability.
 * 3. Invalid-input handling.
 */

#include "../physics/quantum_algorithms.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.10f expected=%.10f err=%.2e\n", label, got, expected,
         err);

  if (err > tol) {
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

static void test_deutsch_jozsa(void) {
  printf("test_deutsch_jozsa:\n");

  dj_result_t r0 = deutsch_jozsa(3, DJ_CONSTANT_0, NULL, 0);
  check_close(r0.p_all_zero, 1.0, 1e-9, "CONSTANT_0: P(all-zero) = 1");
  check_true(r0.is_constant == 1,
             "CONSTANT_0: correctly identified as constant");

  dj_result_t r1 = deutsch_jozsa(3, DJ_CONSTANT_1, NULL, 0);
  check_close(r1.p_all_zero, 1.0, 1e-9, "CONSTANT_1: P(all-zero) = 1");
  check_true(r1.is_constant == 1,
             "CONSTANT_1: correctly identified as constant");

  const int parity1[1] = {0};
  dj_result_t rb1 = deutsch_jozsa(3, DJ_BALANCED, parity1, 1);
  check_close(rb1.p_all_zero, 0.0, 1e-9,
              "BALANCED (single-bit parity): P(all-zero) = 0");
  check_true(rb1.is_constant == 0,
             "BALANCED (single-bit): correctly identified as balanced");

  const int parity2[2] = {0, 1};
  dj_result_t rb2 = deutsch_jozsa(3, DJ_BALANCED, parity2, 2);
  check_close(rb2.p_all_zero, 0.0, 1e-9,
              "BALANCED (2-bit parity): P(all-zero) = 0");
  check_true(rb2.is_constant == 0,
             "BALANCED (2-bit): correctly identified as balanced");
}

static void test_deutsch_jozsa_invalid_input(void) {
  printf("test_deutsch_jozsa_invalid_input:\n");

  dj_result_t r1 = deutsch_jozsa(0, DJ_CONSTANT_0, NULL, 0);
  check_true(r1.p_all_zero < 0.0, "n_input=0 rejected");

  dj_result_t r2 = deutsch_jozsa(3, DJ_BALANCED, NULL, 0);
  check_true(r2.p_all_zero < 0.0, "BALANCED with no parity qubits rejected");
}

static void test_grover(void) {
  printf("test_grover:\n");

  int cases[3][2] = {{3, 5}, {5, 17}, {6, 50}}; // {n_qubits, target}

  for (int c = 0; c < 3; c++) {
    int n_qubits = cases[c][0], target = cases[c][1];
    int N = 1 << n_qubits;

    grover_result_t r = grover_search(n_qubits, target, -1);

    double theta = asin(1.0 / sqrt((double)N));
    double p_theory = pow(sin((2 * r.n_iterations + 1) * theta), 2);

    char label[64];
    snprintf(label, sizeof label, "n=%d target=%d: P(target) matches theory",
             n_qubits, target);
    check_close(r.p_target, p_theory, 1e-6, label);

    double sum_p = 0.0;
    for (int i = 0; i < N; i++) {
      sum_p += r.probabilities[i];
    }
    snprintf(label, sizeof label,
             "n=%d target=%d: probabilities normalize to 1", n_qubits, target);
    check_close(sum_p, 1.0, 1e-9, label);

    free(r.probabilities);
  }
}

static void test_grover_overrotation(void) {
  printf("test_grover_overrotation:\n");

  int n_qubits = 4, target = 7;
  int optimal_k = 3; // round((\pi / 4) * \sqrt(16)) = 3
  int trough_k = 6;  // TEST: verified via the closed-form formula
  // NOTE: Grover's success probability is periodic in k, not monotonic, so
  // "more iterations" isn't worse

  grover_result_t r_optimal = grover_search(n_qubits, target, optimal_k);
  grover_result_t r_trough = grover_search(n_qubits, target, trough_k);

  printf("  optimal k=%d: P(target)=%.4f\n", optimal_k, r_optimal.p_target);
  printf("  trough k=%d: P(target)=%.4f\n", trough_k, r_trough.p_target);

  check_true(r_trough.p_target < 0.1, "a k chosen to sit near a probability "
                                      "trough gives low success probability");
  check_true(r_trough.p_target < r_optimal.p_target,
             "trough k gives lower success probability than optimal k \n"
             "(Grover's success rate is periodic, not monotonically increasing "
             "with more iterations)");

  free(r_optimal.probabilities);
  free(r_trough.probabilities);
}

static void test_grover_invalid_input(void) {
  printf("test_grover_invalid_input:\n");

  grover_result_t r1 = grover_search(0, 0, -1);
  check_true(r1.probabilities == NULL, "n_qubits=0 rejected");

  grover_result_t r2 = grover_search(3, 100, -1); // N=8, target out of range
  check_true(r2.probabilities == NULL, "target out of range rejected");

  grover_result_t r3 = grover_search(3, -1, -1);
  check_true(r3.probabilities == NULL, "negative target rejected");
}

int main(void) {
  test_deutsch_jozsa();
  test_deutsch_jozsa_invalid_input();
  test_grover();
  test_grover_overrotation();
  test_grover_invalid_input();

  if (failures == 0) {
    printf("\nAll test_quantum_algorithms checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
