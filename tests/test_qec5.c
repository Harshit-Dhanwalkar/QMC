/*
 * Test: 5-qubit "perfect" quantum error correction code
 *
 * The code must correctly diagnose (via syndrome lookup, matching the known
 * syndrome table) and exactly correct every possible single-qubit error: no
 * error, or X/Y/Z on any of the 5 data qubits (16 cases total). Unlike the
 * 9-qubit Shor code, the 5-qubit code's syndrome identifies the exact Pauli
 * type directly, so recovery must be exact (no global-phase ambiguity) in every
 * one of the 16 cases.
 */

#include "../core/complex.h"
#include "../physics/qec.h"
#include "../physics/qec5.h"
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

// Independently-derived syndrome table (qubit, Pauli type) for every one of 15
// single-qubit errors, indexed [qubit][type] with type 0=X,1=Y,2=Z;
// {s1,s2,s3,s4} is expected 4-bit syndrome
// NOTE: Cross-checked against physics/qec5.c's qec5_syndrome_table
static const int expected_syndrome[5][3][4] = {
    /* q0 */ {{0, 0, 0, 1}, {1, 0, 1, 1}, {1, 0, 1, 0}},
    /* q1 */ {{1, 0, 0, 0}, {1, 1, 0, 1}, {0, 1, 0, 1}},
    /* q2 */ {{1, 1, 0, 0}, {1, 1, 1, 0}, {0, 0, 1, 0}},
    /* q3 */ {{0, 1, 1, 0}, {1, 1, 1, 1}, {1, 0, 0, 1}},
    /* q4 */ {{0, 0, 1, 1}, {0, 1, 1, 1}, {0, 1, 0, 0}},
};

static const qec_error_type_t types[3] = {QEC_ERROR_X, QEC_ERROR_Y,
                                          QEC_ERROR_Z};
static const char *const type_names[3] = {"X", "Y", "Z"};

static void run_all_errors(complex_t alpha, complex_t beta,
                           const char *state_label) {
  const double u[4] = {0.3, 0.7, 0.2, 0.8};

  // No error
  {
    qec5_result_t r = qec5_run(alpha, beta, -1, QEC_ERROR_X, u);
    char label[96];

    snprintf(label, sizeof label, "%s no error: no correction triggered",
             state_label);
    check_true(r.corrected_qubit == -1 && r.syndrome[0] == 0 &&
                   r.syndrome[1] == 0 && r.syndrome[2] == 0 &&
                   r.syndrome[3] == 0,
               label);

    snprintf(label, sizeof label, "%s no error: recovers exact original state",
             state_label);
    check_state_exact(r.recovered_alpha, r.recovered_beta, alpha, beta, label);
  }

  // Every single-qubit X, Y, and Z error
  for (int q = 0; q < 5; q++) {
    for (int t = 0; t < 3; t++) {
      qec5_result_t r = qec5_run(alpha, beta, q, types[t], u);
      char label[128];

      snprintf(label, sizeof label,
               "%s %s error on q%d: syndrome matches expected table",
               state_label, type_names[t], q);
      check_true(r.syndrome[0] == expected_syndrome[q][t][0] &&
                     r.syndrome[1] == expected_syndrome[q][t][1] &&
                     r.syndrome[2] == expected_syndrome[q][t][2] &&
                     r.syndrome[3] == expected_syndrome[q][t][3],
                 label);

      snprintf(label, sizeof label,
               "%s %s error on q%d: corrected_qubit/type match injected error",
               state_label, type_names[t], q);
      check_true(r.corrected_qubit == q && r.corrected_type == types[t], label);

      snprintf(label, sizeof label,
               "%s %s error on q%d: recovers exact original state (no phase "
               "ambiguity)",
               state_label, type_names[t], q);
      check_state_exact(r.recovered_alpha, r.recovered_beta, alpha, beta,
                        label);
    }
  }
}

int main(void) {
  printf("test_qec5:\n");

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
    printf("\nAll test_qec5 checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
