/*
 * Test: distance-3 rotated surface code ([[9,1,3]]).
 *
 * NOTE: Checks all 28 "no error or single-qubit error" cases (9 qubits * 3
 * Pauli types + no-error) against the syndrome and canonical (qubit, type)
 * correction and verifies exact recovery of the original logical amplitudes in
 * every case - including 4 syndrome-aliased pairs, where two distinct physical
 * errors share a syndrome but the canonical correction still restores the state
 * exactly (they differ by a stabilizer element, not a logical operator).
 */

#include "../core/complex.h"
#include "../physics/qec_surface17.h"
#include "physics/qec.h"
#include <math.h>
#include <stdio.h>

static int failures = 0;

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

static double cabs2_diff(complex_t a, complex_t b) {
  complex_t d = c_sub(a, b);

  return c_abs2(d);
}

static void check_state_exact(complex_t got_a, complex_t got_b, complex_t exp_a,
                              complex_t exp_b, const char *label) {
  double err = sqrt(cabs2_diff(got_a, exp_a) + cabs2_diff(got_b, exp_b));
  printf("  %s: err=%.2e\n", label, err);
  if (err > 1e-9) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static const int expected_syndrome[9][3][8] = {
    /* q0 */
    {{0, 1, 0, 0, 0, 0, 0, 0},
     {0, 1, 1, 0, 0, 0, 0, 0},
     {0, 0, 1, 0, 0, 0, 0, 0}},
    /* q1 */
    {{0, 1, 0, 0, 1, 0, 0, 0},
     {0, 1, 1, 0, 1, 0, 0, 0},
     {0, 0, 1, 0, 0, 0, 0, 0}},
    /* q2 */
    {{0, 0, 0, 0, 1, 0, 0, 0},
     {0, 0, 0, 0, 1, 0, 0, 1},
     {0, 0, 0, 0, 0, 0, 0, 1}},
    /* q3 */
    {{0, 0, 0, 1, 0, 0, 0, 0},
     {1, 0, 1, 1, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0, 0, 0}},
    /* q4 */
    {{0, 0, 0, 1, 1, 0, 0, 0},
     {0, 0, 1, 1, 1, 1, 0, 0},
     {0, 0, 1, 0, 0, 1, 0, 0}},
    /* q5 */
    {{0, 0, 0, 0, 1, 0, 0, 0},
     {0, 0, 0, 0, 1, 1, 0, 1},
     {0, 0, 0, 0, 0, 1, 0, 1}},
    /* q6 */
    {{0, 0, 0, 1, 0, 0, 0, 0},
     {1, 0, 0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0, 0}},
    /* q7 */
    {{0, 0, 0, 1, 0, 0, 1, 0},
     {0, 0, 0, 1, 0, 1, 1, 0},
     {0, 0, 0, 0, 0, 1, 0, 0}},
    /* q8 */
    {{0, 0, 0, 0, 0, 0, 1, 0},
     {0, 0, 0, 0, 0, 1, 1, 0},
     {0, 0, 0, 0, 0, 1, 0, 0}},
};

static const qec_error_type_t types[3] = {QEC_ERROR_X, QEC_ERROR_Y,
                                          QEC_ERROR_Z};
static const char *const type_names[3] = {"X", "Y", "Z"};

static void run_all_errors(complex_t alpha, complex_t beta,
                           const char *state_label) {
  const double u[8] = {0.3, 0.7, 0.2, 0.8, 0.1, 0.9, 0.4, 0.6};

  {
    qec_surface17_result_t r =
        qec_surface17_run(alpha, beta, -1, QEC_ERROR_X, u);
    char label[96];

    snprintf(label, sizeof label, "%s no error: no correction triggered",
             state_label);
    check_true(r.corrected_qubit == -1, label);

    snprintf(label, sizeof label, "%s no error: recovers exact original state",
             state_label);
    check_state_exact(r.recovered_alpha, r.recovered_beta, alpha, beta, label);
  }

  for (int q = 0; q < 9; q++) {
    for (int t = 0; t < 3; t++) {
      qec_surface17_result_t r = qec_surface17_run(alpha, beta, q, types[t], u);
      char label[128];
      int synd_ok = 1;

      for (int k = 0; k < 8; k++) {
        if (r.syndrome[k] != expected_syndrome[q][t][k]) {
          synd_ok = 0;
        }
      }

      snprintf(label, sizeof label,
               "%s %s error on q%d: syndrome matches expected table",
               state_label, type_names[t], q);
      check_true(synd_ok, label);

      snprintf(label, sizeof label,
               "%s %s error on q%d: recovers exact original state", state_label,
               type_names[t], q);
      check_state_exact(r.recovered_alpha, r.recovered_beta, alpha, beta,
                        label);
    }
  }
}

int main(void) {
  printf("test_qec_surface17:\n");

  complex_t alpha1 = c_real(0.6);
  complex_t beta1 = c_new(0.7368487952023082, 0.31153467384692046);
  run_all_errors(alpha1, beta1, "state1");

  complex_t alpha2 = c_new(0.6, 0.2);
  complex_t beta2raw = c_new(0.5, -0.3);
  double norm2 = sqrt(c_abs2(alpha2) + c_abs2(beta2raw));
  complex_t alpha2n = c_scale(alpha2, 1.0 / norm2);
  complex_t beta2n = c_scale(beta2raw, 1.0 / norm2);
  run_all_errors(alpha2n, beta2n, "state2");

  if (failures == 0) {
    printf("\nAll test_qec_surface17 checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
