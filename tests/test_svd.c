/*
 * Test: svd_decompose (core/linalg/svd.c)
 *
 * NOTE: Regression: svd_decompose previously built V's columns via
 * `eig->eigenvectors[j].data[r]`, but eig->eigenvectors is a single n x n
 * cmatrix_t (eigenvectors stored as columns), not an array of per-eigenvector
 * structs
 *
 * Checks (against 4x3 real matrix with singular values):
 *  1. Singular values match the numpy reference.
 *  2. A = U * diag(S) * V^T reconstructs original matrix
 *  3. U and V both have orthonormal columns.
 */

#include "../core/complex.h"
#include "../core/linalg/svd.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  if (err > tol) {
    printf("  FAIL: %s (got=%.10f expected=%.10f err=%.2e > tol=%.2e)\n", label,
           got, expected, err, tol);
    failures++;
  } else {
    printf("  ok:   %s (err=%.2e)\n", label, err);
  }
}

static void test_svd_4x3(void) {
  printf(
      " === svd_decompose: 4x3 real matrix, 3 distinct singular values ===\n");

  int m = 4, n = 3;
  cmatrix_t *A = cmatrix_alloc(m, n);
  static const double vals[4][3] = {
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 10.0}, {1.0, 0.0, 1.0}};

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(A, i, j) = c_real(vals[i][j]);
    }
  }

  cmatrix_t *U = cmatrix_alloc(m, n);
  cmatrix_t *V = cmatrix_alloc(n, n);
  cvector_t *S = cvector_alloc(n);

  int ret = svd_decompose(A, U, S, V);
  check_close(ret, 0.0, 1e-14, "svd_decompose returns success");
  if (ret != 0) {
    return;
  }

  // Reference singular values 
  static const double S_ref[3] = {17.4508955846, 0.9869391657, 0.7015656614};
  for (int i = 0; i < n; i++) {
    char label[64];
    snprintf(label, sizeof label, "singular value S[%d]", i);

    check_close(S->data[i].re, S_ref[i], 1e-6, label);
  }

  // Reconstruction: A =?= U * diag(S) * V^T
  double max_recon_err = 0.0;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      complex_t sum = c_zero();

      for (int k = 0; k < n; k++) {
        complex_t term =
            c_mul(CMAT(U, i, k), c_scale(CMAT(V, j, k), S->data[k].re));

        sum = c_add(sum, term);
      }

      double err = fabs(sum.re - vals[i][j]);
      if (err > max_recon_err) {
        max_recon_err = err;
      }
    }
  }

  check_close(max_recon_err, 0.0, 1e-9,
              "A == U * diag(S) * V^T (reconstruction)");

  // U columns orthonormal: U^T U == I_n
  double max_u_orth_err = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      complex_t dot = c_zero();
      for (int k = 0; k < m; k++) {
        dot = c_add(dot, c_mul(c_conj(CMAT(U, k, i)), CMAT(U, k, j)));
      }

      double expected = (i == j) ? 1.0 : 0.0;
      double err = fabs(dot.re - expected);
      if (err > max_u_orth_err) {
        max_u_orth_err = err;
      }
    }
  }

  check_close(max_u_orth_err, 0.0, 1e-9, "U columns orthonormal");

  // V columns orthonormal: V^T V == I_n
  double max_v_orth_err = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      complex_t dot = c_zero();

      for (int k = 0; k < n; k++) {
        dot = c_add(dot, c_mul(c_conj(CMAT(V, k, i)), CMAT(V, k, j)));
      }

      double expected = (i == j) ? 1.0 : 0.0;
      double err = fabs(dot.re - expected);
      if (err > max_v_orth_err) {
        max_v_orth_err = err;
      }
    }
  }

  check_close(max_v_orth_err, 0.0, 1e-9, "V columns orthonormal");

  cmatrix_free(A);
  cmatrix_free(U);
  cmatrix_free(V);
  cvector_free(S);
}

int main(void) {
  printf(" > SVD decomposition tests\n");

  test_svd_4x3();

  if (failures == 0) {
    printf("\nAll test_svd checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_svd check(s) FAILED.\n", failures);
    return 1;
  }
}
