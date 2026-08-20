/*
 * Test: Quantum Phase Estimation (QPE), physics/qpe.c.
 *
 * Closed-form result (Nielsen & Chuang eq. 5.34): for a single-qubit eigenstate
 * with exact eigenphase \phi, n_count counting qubits, the post-inverse-QFT
 * counting-register distribution is :
 *   P(x) = (1/N^2) * |sum_{k=0}^{N-1} \exp^{2 * \pi * i * k * (\phi - x/N)}|^2,
 *          N=2^n_count
 * which collapses to a delta function (P=1 at x = \phi*N) whenever \phi*N is an
 * integer, and a sinc^2-like spread otherwise.
 *
 * 1. Single-qubit target, phase gate diag(1, \exp^{2 * \pi * i * \phi}) with
 *    eigenstate |1>, phi=3/8 exactly representable with n_count=3 bits: exact
 *    peak, P=1.0.
 * 2. Same setup, phi=0.2 (not exactly representable) with n_count=5 bits:
 *    checks best estimate lands within 1/2^n_count of true phase, and
 *    cross-checks the full distribution against the closed form above.
 * 3. Two-qubit target register, U=diag(1,1,1,e^{2*pi*i*phi}) with eigenstate
 *    |11>, phi=1/4, n_count=2: exercises qstate_apply_controlled_unitary's
 *    multi-qubit-target path (test 1/2 only exercise single-target-qubit case).
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/qpe.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);
  if (err > tol) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAILED");
  if (!cond) {
    failures++;
  }
}

/* Closed-form QPE distribution (Nielsen & Chuang eq. 5.34) for a
   single-eigenstate exact-phase input, computed independently of qpe.c for
   cross-checking the full C-computed distribution, not just its peak. */
static double closed_form_prob(double phi, int n_count, int x) {
  long long N = 1LL << n_count;
  double re = 0.0, im = 0.0;

  for (long long k = 0; k < N; k++) {
    double angle = 2.0 * M_PI * (double)k * (phi - (double)x / (double)N);

    re += cos(angle);
    im += sin(angle);
  }

  return (re * re + im * im) / ((double)N * (double)N);
}

static void test_h2_exact_phase(void) {
  printf("test_qpe_exact_phase_3_8:\n");
  double phi = 3.0 / 8.0;
  int n_count = 3;

  cmatrix_t *U = cmatrix_alloc(2, 2);
  CMAT(U, 0, 0) = c_real(1.0);
  CMAT(U, 0, 1) = c_zero();
  CMAT(U, 1, 0) = c_zero();
  CMAT(U, 1, 1) = c_new(cos(2.0 * M_PI * phi), sin(2.0 * M_PI * phi));

  complex_t eigenstate[2] = {c_zero(), c_real(1.0)}; // |1>

  qpe_distribution_t *d = qpe_estimate_phase(n_count, 1, U, eigenstate);
  check_true(d != NULL, "qpe_estimate_phase allocates");
  if (d) {
    int expected_bin = (int)round(phi * (1 << n_count));

    check_close(d->phi_estimate, phi, 1e-9, "\\phi estimate matches exactly");
    check_close(d->probabilities[expected_bin], 1.0, 1e-9,
                "peak probability is exactly 1 (\\phi exactly representable)");

    qpe_distribution_free(d);
  }

  cmatrix_free(U);
}

static void test_inexact_phase_vs_closed_form(void) {
  printf("test_qpe_inexact_phase_vs_closed_form:\n");
  double phi = 0.2;
  int n_count = 5;

  cmatrix_t *U = cmatrix_alloc(2, 2);
  CMAT(U, 0, 0) = c_real(1.0);
  CMAT(U, 0, 1) = c_zero();
  CMAT(U, 1, 0) = c_zero();
  CMAT(U, 1, 1) = c_new(cos(2.0 * M_PI * phi), sin(2.0 * M_PI * phi));

  complex_t eigenstate[2] = {c_zero(), c_real(1.0)};

  qpe_distribution_t *d = qpe_estimate_phase(n_count, 1, U, eigenstate);
  check_true(d != NULL, "qpe_estimate_phase allocates");
  if (d) {
    check_close(d->phi_estimate, phi, 1.0 / (1 << n_count),
                "\\phi estimate within 1/2^n_count of true phase");

    double max_dist_err = 0.0;
    long long n_bins = 1LL << n_count;

    for (long long x = 0; x < n_bins; x++) {
      double expected = closed_form_prob(phi, n_count, (int)x);
      double err = fabs(d->probabilities[x] - expected);

      if (err > max_dist_err) {
        max_dist_err = err;
      }
    }

    check_close(max_dist_err, 0.0, 1e-9,
                "full distribution matches closed form at every bin");

    double total = 0.0;
    for (long long x = 0; x < n_bins; x++) {
      total += d->probabilities[x];
    }

    check_close(total, 1.0, 1e-9, "probabilities sum to 1");

    qpe_distribution_free(d);
  }

  cmatrix_free(U);
}

static void test_two_qubit_target_register(void) {
  printf("test_qpe_two_qubit_target:\n");
  /* U = diag(1, 1, 1, \exp^{2 * \pi * i * \phi}) on 2 qubits, eigenstate |11>:
   * exercises multi-qubit-target path of qstate_apply_controlled_unitary. */
  double phi = 0.25;
  int n_count = 2;

  cmatrix_t *U = cmatrix_alloc(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      CMAT(U, i, j) = c_zero();
    }
  }

  CMAT(U, 0, 0) = c_real(1.0);
  CMAT(U, 1, 1) = c_real(1.0);
  CMAT(U, 2, 2) = c_real(1.0);
  CMAT(U, 3, 3) = c_new(cos(2.0 * M_PI * phi), sin(2.0 * M_PI * phi));

  complex_t eigenstate[4] = {c_zero(), c_zero(), c_zero(), c_real(1.0)}; // |11>

  qpe_distribution_t *d = qpe_estimate_phase(n_count, 2, U, eigenstate);
  check_true(d != NULL, "qpe_estimate_phase allocates (2-qubit target)");
  if (d) {
    int expected_bin = (int)round(phi * (1 << n_count));

    check_close(d->phi_estimate, phi, 1e-9,
                "phi estimate matches exactly (2-qubit target)");
    check_close(d->probabilities[expected_bin], 1.0, 1e-9,
                "peak probability is exactly 1 (2-qubit target)");

    qpe_distribution_free(d);
  }

  cmatrix_free(U);
}

int main(void) {
  test_h2_exact_phase();
  test_inexact_phase_vs_closed_form();
  test_two_qubit_target_register();

  if (failures > 0) {
    printf("\n%d test_qpe check(s) FAILED.\n", failures);
    return 1;
  }
  printf("\nAll test_qpe checks passed.\n");
  return 0;
}
