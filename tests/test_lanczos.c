/*
Test: Lanczos iteration for lowest eigenvalues of a Hermitian sparse matrix.

1. Diagonal matrix sanity check: for a diagonal Hermitian matrix, k lowest
   eigenvalues are trivially k smallest diagonal entries, and eigenvectors are
   corresponding coordinate vectors. Deterministic, exact (to numerical
   tolerance).
2. Cross-check against the project's already-trusted dense Hermitian eigensolver
   (cmatrix_eigh): build a random Hermitian sparse matrix, run both solvers,
   compare k lowest eigenvalues and verify Ax = \lambda * x for Lanczos
   eigenvectors directly.
3. Eigenvector orthonormality: <v_i, v_j> = \delta_ij for returned eigenvectors.
4. Invalid-input handling.
*/

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../core/random.h"
#include "../core/sparse.h"
#include "../core/vector.h"
#include <math.h>
#include <stdint.h>
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
  printf("test_diagonal_matrix:\n");

  int n = 8;
  double diag_vals[8] = {5.0, 1.0, 8.0, 0.5, 3.0, 9.0, 2.0, 6.0};

  cmatrix_t *dense = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(dense, i, j) = c_zero();
    }

    CMAT(dense, i, i) = c_real(diag_vals[i]);
  }

  sparse_matrix_t *A = sparse_from_dense(dense, 1e-12);
  cmatrix_free(dense);

  int k = 3;
  lanczos_result_t *res = lanczos_eigs(A, k, 30, 1e-10);

  check_true(res != NULL, "lanczos_eigs succeeds on diagonal matrix");
  if (res) {
    // sorted smallest diagonal entries: 0.5, 1.0, 2.0
    const double expected[3] = {0.5, 1.0, 2.0};
    for (int i = 0; i < k; i++) {
      char label[32];
      snprintf(label, sizeof label, "eigenvalue[%d]", i);

      check_close(res->values[i], expected[i], 1e-8, label);
    }

    lanczos_free(res);
  }

  sparse_free(A);
}

// Builds a random Hermitian dense matrix (small, size n) with fixed seed
static cmatrix_t *random_hermitian(int n, uint64_t seed) {
  rng_state_t rng;
  rng_seed(&rng, seed);

  cmatrix_t *A = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    double re = rng_uniform_range(&rng, -2.0, 2.0);
    CMAT(A, i, i) = c_real(re);

    for (int j = i + 1; j < n; j++) {
      double re_ij = rng_uniform_range(&rng, -1.0, 1.0);
      double im_ij = rng_uniform_range(&rng, -1.0, 1.0);
      CMAT(A, i, j) = c_new(re_ij, im_ij);
      CMAT(A, j, i) = c_conj(CMAT(A, i, j));
    }
  }

  return A;
}

static void test_random_hermitian_vs_dense(void) {
  printf("test_random_hermitian_vs_dense:\n");

  int n = 12;
  cmatrix_t *dense = random_hermitian(n, 20260802ULL);
  cmatrix_t *dense_for_ref = cmatrix_copy(dense);

  eigen_t *dense_eig = cmatrix_eigh_complex(dense_for_ref);
  cmatrix_free(dense_for_ref);

  sparse_matrix_t *A = sparse_from_dense(dense, 1e-12);

  int k = 3;
  lanczos_result_t *res = lanczos_eigs(A, k, 40, 1e-10);

  check_true(res != NULL && dense_eig != NULL,
             "both solvers succeed on random Hermitian matrix");

  if (res && dense_eig) {
    for (int i = 0; i < k; i++) {
      char label[48];
      snprintf(label, sizeof label, "eigenvalue[%d] vs dense cmatrix_eigh", i);
      check_close(res->values[i], dense_eig->eigenvalues[i], 1e-6, label);
    }

    // Direct residual check: A * v_i - \lambda_i * v_i should be ~0
    cvector_t *v = cvector_alloc(n);
    cvector_t *Av = cvector_alloc(n);
    for (int i = 0; i < k; i++) {
      for (int row = 0; row < n; row++) {
        v->data[row] = CMAT(res->vectors, row, i);
      }
      sparse_mv(A, v, Av);

      double residual = 0.0;
      for (int row = 0; row < n; row++) {
        complex_t diff =
            c_sub(Av->data[row], c_scale(v->data[row], res->values[i]));
        residual += c_abs2(diff);
      }
      residual = sqrt(residual);

      char label[48];
      snprintf(label, sizeof label, "||A*v[%d] - lambda[%d]*v[%d]|| ~ 0", i, i,
               i);
      check_close(residual, 0.0, 1e-6, label);
    }
    cvector_free(v);
    cvector_free(Av);

    // Orthonormality of returned eigenvectors
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        complex_t dot = c_zero();

        for (int row = 0; row < n; row++) {
          dot = c_add(dot, c_mul(c_conj(CMAT(res->vectors, row, i)),
                                 CMAT(res->vectors, row, j)));
        }

        const double expected_re = (i == j) ? 1.0 : 0.0;
        char label[48];
        snprintf(label, sizeof label, "<v[%d],v[%d]>.re", i, j);
        check_close(dot.re, expected_re, 1e-6, label);

        char label2[48];
        snprintf(label2, sizeof label2, "<v[%d],v[%d]>.im", i, j);
        check_close(dot.im, 0.0, 1e-6, label2);
      }
    }
  }

  if (res) {
    lanczos_free(res);
  }
  if (dense_eig) {
    eigen_free(dense_eig);
  }
  cmatrix_free(dense);
  sparse_free(A);
}

static void test_invalid_input(void) {
  printf("test_invalid_input:\n");

  cmatrix_t *dense = cmatrix_alloc(3, 3);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      CMAT(dense, i, j) = c_zero();
    }

    CMAT(dense, i, i) = c_real(1.0 + i);
  }

  sparse_matrix_t *A = sparse_from_dense(dense, 1e-12);
  cmatrix_free(dense);

  check_true(lanczos_eigs(NULL, 1, 10, 1e-10) == NULL, "NULL matrix rejected");
  check_true(lanczos_eigs(A, 0, 10, 1e-10) == NULL, "k=0 rejected");
  check_true(lanczos_eigs(A, 10, 10, 1e-10) == NULL, "k > n rejected (n=3, k=10)");
  check_true(lanczos_eigs(A, 3, 1, 1e-10) == NULL, "max_iter < k rejected");
  check_true(lanczos_eigs(A, 1, 10, 0.0) == NULL, "tol<=0 rejected");
  check_true(lanczos_eigs(A, 1, 10, -1.0) == NULL, "negative tol rejected");

  sparse_free(A);
}

int main(void) {
  test_diagonal_matrix();
  test_random_hermitian_vs_dense();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_lanczos checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
