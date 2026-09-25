/*
 * Test: shift-invert eigensolver (inverse iteration on A - \sigma * I via
 * a single dense LU factorization, with Gram-Schmidt deflation for k > 1)
 *
 * 1. Diagonal matrix sanity check: for a diagonal Hermitian matrix,
 *    eigenvalues nearest a chosen shift are trivially identifiable by
 *    inspection. Deterministic, exact to numerical tolerance
 * 2. Cross-check against dense Hermitian  eigensolver (cmatrix_eigh): build a
 *    random Hermitian sparse matrix, shift near an interior eigenvalue, and
 *    verify shift_invert_eigs recovers same eigenvalues dense solver reports as
 *    nearest to \sigma, plus A v = \lambda v directly for each returned
 *    eigenvector
 * 3. Eigenvector orthonormality across k returned eigenvectors
 * 4. Invalid-input handling, including exact-eigenvalue (singular  shift) case
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../core/random.h"
#include "../core/shift_invert.h"
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

static void test_diagonal_matrix(void) {
  printf("  === Test diagonal matrix ===\n");

  /* Diagonal entries 0.5, 1.0, 2.0, 5.0, 9.0. Shift near 2.0: two nearest
   * eigenvalues are 1.0 and 2.0 */
  int n = 5;
  const double diag[5] = {9.0, 0.5, 5.0, 2.0, 1.0};
  cmatrix_t *dense = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    CMAT(dense, i, i) = c_real(diag[i]);
  }

  sparse_matrix_t *A = sparse_from_dense(dense, 1e-12);
  cmatrix_free(dense);

  shift_invert_result_t *res = shift_invert_eigs(A, 2.1, 2, 50, 1e-12);
  check_true(res != NULL, "shift_invert_eigs succeeds on diagonal matrix");
  if (res) {
    check_true(res->n == 2, "res->n equals k, not matrix dimension");
    check_close(res->values[0], 1.0, 1e-8, "nearest eigenvalue to \\sigma=2.1");
    check_close(res->values[1], 2.0, 1e-8, "2nd-nearest eigenvalue");
  }

  shift_invert_result_free(res);
  sparse_free(A);
}

static void test_random_hermitian_vs_dense(void) {
  printf("  === Test Random Hermitian vs Dense ===\n");

  int n = 10;
  rng_state_t rng;
  rng_seed(&rng, 0x517E1F0ULL);

  cmatrix_t *H = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      complex_t z =
          c_new(rng_gaussian(&rng), (i == j) ? 0.0 : rng_gaussian(&rng));

      CMAT(H, i, j) = z;
      CMAT(H, j, i) = c_conj(z);
    }
  }

  sparse_matrix_t *A = sparse_from_dense(H, 1e-12);

  cmatrix_t *H_for_ref = cmatrix_copy(H);
  eigen_t *dense_eig = cmatrix_eigh_complex(H_for_ref);
  cmatrix_free(H_for_ref);
  check_true(dense_eig != NULL, "reference dense eigensolver succeeds");

  if (dense_eig) {
    const double *dense_evals = dense_eig->eigenvalues;

    // Shift near middle eigenvalue; ask for 2 nearest
    double sigma = dense_evals[n / 2] + 0.05;

    shift_invert_result_t *res = shift_invert_eigs(A, sigma, 2, 200, 1e-12);
    check_true(res != NULL, "shift_invert_eigs succeeds");

    if (res) {
      // Find 2 dense eigenvalues actually nearest \sigma, for comparison
      int best[2] = {-1, -1};

      for (int pick = 0; pick < 2; pick++) {
        double best_dist = 1e300;

        for (int i = 0; i < n; i++) {
          if (i == best[0]) {
            continue;
          }

          double dist = fabs(dense_evals[i] - sigma);
          if (dist < best_dist) {
            best_dist = dist;
            best[pick] = i;
          }
        }
      }

      if (best[0] > best[1]) {
        int tmp = best[0];
        best[0] = best[1];
        best[1] = tmp;
      }

      check_close(res->values[0], dense_evals[best[0]], 1e-6,
                  "1st shift-invert eigenvalue matches dense reference");
      check_close(res->values[1], dense_evals[best[1]], 1e-6,
                  "2nd shift-invert eigenvalue matches dense reference");

      // A v = \lambda v directly, for each returned eigenvector
      cvector_t *v = cvector_alloc(n);
      cvector_t *Av = cvector_alloc(n);
      for (int col = 0; col < res->n; col++) {
        for (int i = 0; i < n; i++) {
          v->data[i] = CMAT(res->vectors, i, col);
        }

        sparse_mv(A, v, Av);

        double max_resid = 0.0;
        for (int i = 0; i < n; i++) {
          complex_t expected = c_scale(v->data[i], res->values[col]);

          double err = c_abs(c_sub(Av->data[i], expected));
          if (err > max_resid) {
            max_resid = err;
          }
        }

        char label[64];
        snprintf(label, sizeof label, "A v = \\lambda v residual (eigvec %d)",
                 col);
        check_close(max_resid, 0.0, 1e-6, label);
      }

      // Orthonormality across returned eigenvectors
      for (int a = 0; a < res->n; a++) {
        for (int i = 0; i < n; i++) {
          v->data[i] = CMAT(res->vectors, i, a);
        }

        for (int b = 0; b < res->n; b++) {
          cvector_t *w = cvector_alloc(n);

          for (int i = 0; i < n; i++) {
            w->data[i] = CMAT(res->vectors, i, b);
          }

          complex_t ip = cvector_dot(v, w);
          double expected = (a == b) ? 1.0 : 0.0;
          char label[64];
          snprintf(label, sizeof label, "<v%d|v%d> orthonormality", a, b);
          check_close(c_abs(ip), expected, 1e-6, label);
          cvector_free(w);
        }
      }

      cvector_free(v);
      cvector_free(Av);
    }

    shift_invert_result_free(res);
  }

  eigen_free(dense_eig);
  cmatrix_free(H);
  sparse_free(A);
}

static void test_invalid_input(void) {
  printf("  === Test invalid input ===\n");

  int n = 3;
  cmatrix_t *dense = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    CMAT(dense, i, i) = c_real(1.0 + i);
  }

  sparse_matrix_t *A = sparse_from_dense(dense, 1e-12);
  cmatrix_free(dense);

  check_true(shift_invert_eigs(NULL, 0.5, 1, 20, 1e-10) == NULL,
             "NULL matrix rejected");
  check_true(shift_invert_eigs(A, 0.5, 0, 20, 1e-10) == NULL, "k=0 rejected");
  check_true(shift_invert_eigs(A, 0.5, 10, 20, 1e-10) == NULL,
             "k > n rejected (n=3, k=10)");
  check_true(shift_invert_eigs(A, 0.5, 1, 0, 1e-10) == NULL,
             "max_iter=0 rejected");
  check_true(shift_invert_eigs(A, 0.5, 1, 20, 0.0) == NULL, "tol<=0 rejected");
  check_true(shift_invert_eigs(A, 0.5, 1, 20, -1.0) == NULL,
             "negative tol rejected");

  // \\sigma exactly on an eigenvalue: (A - \sigma * I) is exactly singular
  shift_invert_result_t *singular = shift_invert_eigs(A, 2.0, 1, 20, 1e-10);
  check_true(singular == NULL, "exact eigenvalue shift rejected (singular)");
  shift_invert_result_free(singular);

  sparse_free(A);
}

int main(void) {
  test_diagonal_matrix();
  test_random_hermitian_vs_dense();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_shift_invert checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
