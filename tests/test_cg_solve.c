/*
 * Test: Conjugate Gradient solver for Hermitian positive-definite sparse linear
 * systems (A x = b)
 *
 * 1. Discrete-Laplacian sanity check: a tridiagonal (2, -1, -1) matrix is a
 *    theoretical Hermitian positive-definite system; recovers a known exact
 *    solution to near machine precision
 * 2. Cross-check against dense LU solver (lu_decompose/lu_solve): build a
 *    random Hermitian positive-definite sparse matrix (A = B^H B + eps*I),
 * solve same system both ways,  compare solutions and verify residual directly
 * 3. Invalid-input handling
 */

#include "../core/complex.h"
#include "../core/linalg/lu.h"
#include "../core/matrix.h"
#include "../core/random.h"
#include "../core/sparse.h"
#include "../core/vector.h"
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

static void test_discrete_laplacian(void) {
  printf("  === Test discrete laplacian ===\n");

  int n = 10;
  cmatrix_t *dense = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    CMAT(dense, i, i) = c_real(2.0);

    if (i > 0) {
      CMAT(dense, i, i - 1) = c_real(-1.0);
    }
    if (i < n - 1) {
      CMAT(dense, i, i + 1) = c_real(-1.0);
    }
  }

  sparse_matrix_t *A = sparse_from_dense(dense, 1e-12);
  cmatrix_free(dense);

  rng_state_t rng;
  rng_seed(&rng, 0xC65015EULL);

  cvector_t *x_true = cvector_alloc(n);
  for (int i = 0; i < n; i++) {
    x_true->data[i] = c_real(rng_gaussian(&rng));
  }

  cvector_t *b = cvector_alloc(n);
  sparse_mv(A, x_true, b);

  cvector_t *x = cvector_alloc(n);
  cvector_fill(x, c_zero()); // zero initial guess

  int status = cg_solve(A, b, x, n, 1e-14);
  check_true(status == 0, "cg_solve converges");

  double max_err = 0.0;
  for (int i = 0; i < n; i++) {
    double err = c_abs(c_sub(x->data[i], x_true->data[i]));
    if (err > max_err) {
      max_err = err;
    }
  }

  check_close(max_err, 0.0, 1e-9, "max component error vs known solution");

  cvector_free(x_true);
  cvector_free(b);
  cvector_free(x);
  sparse_free(A);
}

static void test_random_hermitian_pd_vs_lu(void) {
  printf("  === Test Random hermitian PD vs LU ===\n");

  int n = 12;
  rng_state_t rng;
  rng_seed(&rng, 0x5EED5EEDULL);

  cmatrix_t *B = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(B, i, j) = c_new(rng_gaussian(&rng), rng_gaussian(&rng));
    }
  }

  cmatrix_t *Badj = cmatrix_adjoint(B);
  cmatrix_t *A_dense = cmatrix_multiply(Badj, B); // B^H B: Hermitian PSD
  for (int i = 0; i < n; i++) {
    CMAT(A_dense, i, i) = c_add(CMAT(A_dense, i, i), c_real(1.0)); // +I: PD
  }

  cmatrix_free(B);
  cmatrix_free(Badj);

  sparse_matrix_t *A = sparse_from_dense(A_dense, 1e-12);

  cvector_t *b = cvector_alloc(n);
  for (int i = 0; i < n; i++) {
    b->data[i] = c_new(rng_gaussian(&rng), rng_gaussian(&rng));
  }

  // Trusted reference: dense LU solve on an independent copy of A_dense
  cmatrix_t *A_lu = cmatrix_copy(A_dense);
  int *pivot = lu_decompose(A_lu);
  cvector_t *x_lu = cvector_alloc(n);
  check_true(pivot != NULL && lu_solve(A_lu, pivot, b, x_lu) == 0,
             "reference LU solve succeeds");

  cvector_t *x_cg = cvector_alloc(n);
  cvector_fill(x_cg, c_zero());
  int status = cg_solve(A, b, x_cg, 200, 1e-12);
  check_true(status == 0, "cg_solve converges");

  double max_err = 0.0;
  for (int i = 0; i < n; i++) {
    double err = c_abs(c_sub(x_cg->data[i], x_lu->data[i]));
    if (err > max_err) {
      max_err = err;
    }
  }

  check_close(max_err, 0.0, 1e-7, "max component error vs dense LU solve");

  // Direct residual check, independent of either solver
  cvector_t *Ax = cvector_alloc(n);
  sparse_mv(A, x_cg, Ax);
  double residual = 0.0;

  for (int i = 0; i < n; i++) {
    double err = c_abs(c_sub(Ax->data[i], b->data[i]));
    if (err > residual) {
      residual = err;
    }
  }

  check_close(residual, 0.0, 1e-7, "max |Ax - b| residual");

  cmatrix_free(A_dense);
  cmatrix_free(A_lu);
  free(pivot);
  sparse_free(A);
  cvector_free(b);
  cvector_free(x_lu);
  cvector_free(x_cg);
  cvector_free(Ax);
}

static void test_invalid_input(void) {
  printf("  === Test invalid input ===\n");

  cmatrix_t *dense = cmatrix_alloc(3, 3);
  for (int i = 0; i < 3; i++) {
    CMAT(dense, i, i) = c_real(1.0 + i);
  }
  sparse_matrix_t *A = sparse_from_dense(dense, 1e-12);
  cmatrix_free(dense);

  cvector_t *b = cvector_alloc(3);
  cvector_t *x = cvector_alloc(3);
  cvector_t *b_wrong = cvector_alloc(4);
  cvector_fill(b, c_one());
  cvector_fill(x, c_zero());

  check_true(cg_solve(NULL, b, x, 10, 1e-10) == -1, "NULL matrix rejected");
  check_true(cg_solve(A, NULL, x, 10, 1e-10) == -1, "NULL b rejected");
  check_true(cg_solve(A, b, NULL, 10, 1e-10) == -1, "NULL x rejected");
  check_true(cg_solve(A, b_wrong, x, 10, 1e-10) == -1,
             "mismatched b size rejected");
  check_true(cg_solve(A, b, x, 0, 1e-10) == -1, "max_iter=0 rejected");
  check_true(cg_solve(A, b, x, 10, 0.0) == -1, "tol<=0 rejected");
  check_true(cg_solve(A, b, x, 10, -1.0) == -1, "negative tol rejected");

  sparse_free(A);
  cvector_free(b);
  cvector_free(x);
  cvector_free(b_wrong);
}

int main(void) {
  test_discrete_laplacian();
  test_random_hermitian_pd_vs_lu();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_cg_solve checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
