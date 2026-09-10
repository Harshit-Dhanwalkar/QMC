/*
 * Test: 3-qubit bit-flip and phase-flip quantum error correction codes.
 *
 * For both codes, every possible single-qubit error (none, or on data qubit 0,
 * 1, or 2) must be correctly diagnosed by syndrome (matching the known syndrome
 * table) and exactly corrected : then recovered logical qubit must equal to
 * original (\alpha, \beta) to machine precision, regardless of which qubit (if
 * any) was hit.
 */

#include "../core/complex.h"
#include "../physics/qec.h"
#include <math.h>
#include <stdio.h>

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

static void run_code_all_errors(qec_code_t code, const char *name) {
  printf("test_qec_%s:\n", name);

  complex_t alpha = c_real(0.6);
  complex_t beta = c_new(0.7368487952023082, 0.31153467384692046);

  const int expected_syndrome[4][2] = {
      {0, 0}, // no error
      {1, 0}, // qubit 0
      {1, 1}, // qubit 1
      {0, 1}, // qubit 2
  };

  for (int e = -1; e < 3; e++) {
    qec_result_t r = qec_run(code, alpha, beta, e, 0.3, 0.7);

    int idx = e + 1;
    char label_syn[64], label_corr[64], label_a[64], label_b[64];

    snprintf(label_syn, sizeof label_syn,
             "error_qubit=%d: syndrome=(%d,%d) matches expected", e,
             r.syndrome_s3, r.syndrome_s4);
    check_true(r.syndrome_s3 == expected_syndrome[idx][0] &&
                   r.syndrome_s4 == expected_syndrome[idx][1],
               label_syn);

    snprintf(label_corr, sizeof label_corr,
             "error_qubit=%d: corrected_qubit=%d matches injected qubit", e,
             r.corrected_qubit);
    check_true(r.corrected_qubit == e, label_corr);

    snprintf(label_a, sizeof label_a, "error_qubit=%d: recovered \\alpha.re",
             e);
    snprintf(label_b, sizeof label_b, "error_qubit=%d: recovered \\beta", e);
    check_close(r.recovered_alpha.re, alpha.re, 1e-9, label_a);
    check_close(r.recovered_beta.re, beta.re, 1e-9, label_b);
    check_close(r.recovered_beta.im, beta.im, 1e-9, label_b);
  }
}

// 9-qubit Shor code: corrects an arbitrary single-qubit error (X, Y, or Z on
// any of the 9 physical qubits, plus the no-error case)
//
// WARN: Y errors: correcting a Y error means the bit-flip syndrome fires (as if
// X) AND the phase-flip syndrome fires (as if Z) independently, so the applied
// correction is Z*X, not Y. Since Y=iXZ, Z*X*Y = i*Identity - the recovered
// state is correct up to an unavoidable, physically unobservable global phase
// of i. This is expected, standard behavior for independent Pauli-frame
// correction, so Y-error cases are checked up to global phase; X, Z, and
// no-error cases recover the exact original alpha/beta with zero phase
// ambiguity and are checked directly.
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

static void check_state_up_to_phase(complex_t got_a, complex_t got_b,
                                    complex_t exp_a, complex_t exp_b,
                                    const char *label) {
  // Divide out whichever expected component has larger magnitude to find global
  // phase, then check both components against that phase times the expected
  // amplitudes
  int use_a = c_abs2(exp_a) >= c_abs2(exp_b);
  complex_t got_ref = use_a ? got_a : got_b;
  complex_t exp_ref = use_a ? exp_a : exp_b;
  double denom = c_abs2(exp_ref);
  complex_t phase = c_scale(c_mul(got_ref, c_conj(exp_ref)), 1.0 / denom);

  double phase_mod_err = fabs(sqrt(c_abs2(phase)) - 1.0);
  complex_t pred_a = c_mul(phase, exp_a);
  complex_t pred_b = c_mul(phase, exp_b);
  double err = sqrt(cabs2_diff(got_a, pred_a) + cabs2_diff(got_b, pred_b));

  printf("  %s: |phase|-1=%.2e state_err=%.2e\n", label, phase_mod_err, err);
  if (phase_mod_err > 1e-6 || err > 1e-8) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void test_shor_code(void) {
  printf("test_qec_shor:\n");

  complex_t alpha = c_new(0.6, 0.2);
  complex_t beta = c_new(0.5, -0.3);
  double norm = sqrt(c_abs2(alpha) + c_abs2(beta));
  alpha = c_scale(alpha, 1.0 / norm);
  beta = c_scale(beta, 1.0 / norm);

  const double u[8] = {0.3, 0.7, 0.2, 0.8, 0.4, 0.6, 0.35, 0.65};

  // No error
  {
    qec_shor_result_t r = qec_shor_run(alpha, beta, -1, QEC_ERROR_X, u);
    check_true(r.block_corrected[0] == -1 && r.block_corrected[1] == -1 &&
                   r.block_corrected[2] == -1 && r.phase_block_corrected == -1,
               "no error: no correction triggered");
    check_state_exact(r.recovered_alpha, r.recovered_beta, alpha, beta,
                      "no error: recovers exact original state");
  }

  // Every single-qubit X error (bit-flip only - phase syndrome must stay
  // silent)
  for (int q = 0; q < 9; q++) {
    qec_shor_result_t r = qec_shor_run(alpha, beta, q, QEC_ERROR_X, u);
    char label[96];

    snprintf(label, sizeof label, "X error on qubit %d: recovers exact state",
             q);
    check_state_exact(r.recovered_alpha, r.recovered_beta, alpha, beta, label);

    snprintf(label, sizeof label,
             "X error on qubit %d: phase syndrome stayed silent", q);
    check_true(r.phase_block_corrected == -1, label);
  }

  // Every single-qubit Z error (phase-flip only - bit-flip syndromes must stay
  // silent)
  for (int q = 0; q < 9; q++) {
    qec_shor_result_t r = qec_shor_run(alpha, beta, q, QEC_ERROR_Z, u);
    char label[96];

    snprintf(label, sizeof label, "Z error on qubit %d: recovers exact state",
             q);
    check_state_exact(r.recovered_alpha, r.recovered_beta, alpha, beta, label);

    snprintf(label, sizeof label,
             "Z error on qubit %d: bit-flip syndromes stayed silent", q);
    check_true(r.block_corrected[0] == -1 && r.block_corrected[1] == -1 &&
                   r.block_corrected[2] == -1,
               label);
  }

  // Every single-qubit Y error: both syndromes fire (bit-flip and phase-flip),
  // recovering the state up to the expected global phase of i.
  for (int q = 0; q < 9; q++) {
    qec_shor_result_t r = qec_shor_run(alpha, beta, q, QEC_ERROR_Y, u);
    char label[96];

    snprintf(label, sizeof label,
             "Y error on qubit %d: recovers state up to global phase", q);
    check_state_up_to_phase(r.recovered_alpha, r.recovered_beta, alpha, beta,
                            label);

    snprintf(label, sizeof label,
             "Y error on qubit %d: both bit-flip and phase syndromes fired", q);
    int block = q / 3;
    check_true(r.block_corrected[block] == q &&
                   r.phase_block_corrected == block,
               label);
  }
}

int main(void) {
  run_code_all_errors(QEC_BITFLIP, "bitflip");
  run_code_all_errors(QEC_PHASEFLIP, "phaseflip");
  test_shor_code();

  if (failures == 0) {
    printf("\nAll test_qec checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
