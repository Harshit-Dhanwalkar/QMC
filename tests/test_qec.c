/*
Test: 3-qubit bit-flip and phase-flip quantum error correction codes.

For both codes, every possible single-qubit error (none, or on data qubit 0, 1,
or 2) must be correctly diagnosed by syndrome (matching the known syndrome
table) and exactly corrected : then the recovered logical qubit must equal to
the original (\alpha, \beta) to machine precision, regardless of which qubit (if
any) was hit.
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

    snprintf(label_a, sizeof label_a, "error_qubit=%d: recovered \\alpha.re", e);
    snprintf(label_b, sizeof label_b, "error_qubit=%d: recovered \\beta", e);
    check_close(r.recovered_alpha.re, alpha.re, 1e-9, label_a);
    check_close(r.recovered_beta.re, beta.re, 1e-9, label_b);
    check_close(r.recovered_beta.im, beta.im, 1e-9, label_b);
  }
}

int main(void) {
  run_code_all_errors(QEC_BITFLIP, "bitflip");
  run_code_all_errors(QEC_PHASEFLIP, "phaseflip");

  if (failures == 0) {
    printf("\nAll test_qec checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
