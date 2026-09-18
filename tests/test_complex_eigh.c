/*
 * Test: complex-Hermitian eigensolver (cmatrix_eigh_complex)
 *
 * 1. 2x2: Pauli \sigma_y matrix has known eigenvalues +-1 and known
 *   eigenvectors (1,+-i) / \sqrt2
 * 2. General NxN Hermitian check:
 *   verify ||H v_k - \lambda_k v_k|| is small for every eigenpair, eigenvectors
 *   are orthonormal (V^\dagger V ~= I). This checks defining property of an
 *   eigendecomposition directly rather than needing a closed-form reference for
 *   arbitrary matrix
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/angular.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check(int cond, const char *msg) {
  if (!cond) {
    printf("  FAIL: %s\n", msg);
    failures++;
  }
}

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  if (fabs(got - expected) > tol) {
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

static void test_sigma_y(void) {
  printf("  === Test: Pauli sigma_y (known eigenvalues +-1) ===\n");

  cmatrix_t *H = cmatrix_alloc(2, 2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      CMAT(H, i, j) = sigma_y[i * 2 + j];
    }
  }

  eigen_t *eig = cmatrix_eigh_complex(H);
  cmatrix_free(H);
  if (!eig) {
    check(0, "cmatrix_eigh_complex returned NULL");
    return;
  }

  const double expected[2] = {-1.0, 1.0};
  for (int k = 0; k < 2; k++) {
    char label[64];
    snprintf(label, sizeof label, "eigenvalue[%d] matches exact", k);
    check_close(eig->eigenvalues[k], expected[k], 1e-8, label);
  }

  // Eigenvector check up to phase: |v . conj(v_expected)| should be ~1.
  // v_expected for \lambda=-1 is (1,-i)/\sqrt(2), for \lambda=+1 is (1,i) /
  // \sqrt(2)
  double inv_sqrt2 = 1.0 / sqrt(2.0);
  complex_t v_exp[2][2] = {{c_real(inv_sqrt2), c_imag(-inv_sqrt2)},
                           {c_real(inv_sqrt2), c_imag(inv_sqrt2)}};

  for (int k = 0; k < 2; k++) {
    complex_t overlap = c_zero();
    for (int i = 0; i < 2; i++) {
      complex_t vi = CMAT(eig->eigenvectors, i, k);
      overlap = c_add(overlap, c_mul(c_conj(vi), v_exp[k][i]));
    }

    char label[64];
    snprintf(label, sizeof label, "eigenvector[%d] overlap magnitude ~ 1", k);
    check_close(sqrt(c_abs2(overlap)), 1.0, 1e-6, label);
  }

  eigen_free(eig);
}

static void test_random_hermitian(void) {
  printf("  === Test: random 5x5 Hermitian (residual + orthonormality) ===\n");

  int n = 5;
  srand(42);
  cmatrix_t *M = cmatrix_alloc(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double re = (double)rand() / RAND_MAX - 0.5;
      double im = (double)rand() / RAND_MAX - 0.5;

      CMAT(M, i, j) = c_add(c_real(re), c_imag(im));
    }
  }

  // Hermitian-symmetrize: H = (M + M^\dagger)/2
  cmatrix_t *H = cmatrix_alloc(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      complex_t sum = c_add(CMAT(M, i, j), c_conj(CMAT(M, j, i)));

      CMAT(H, i, j) = c_scale(sum, 0.5);
    }
  }

  cmatrix_free(M);

  eigen_t *eig = cmatrix_eigh_complex(H);
  if (!eig) {
    check(0, "cmatrix_eigh_complex returned NULL");
    cmatrix_free(H);

    return;
  }

  double tol = 1e-6;

  // Residual check: ||H v_k - \lambda_k v_k|| for each k
  for (int k = 0; k < n; k++) {
    double lambda = eig->eigenvalues[k];
    double residual = 0.0;

    for (int i = 0; i < n; i++) {
      complex_t Hv_i = c_zero();

      for (int j = 0; j < n; j++) {
        Hv_i = c_add(Hv_i, c_mul(CMAT(H, i, j), CMAT(eig->eigenvectors, j, k)));
      }

      complex_t lv_i = c_scale(CMAT(eig->eigenvectors, i, k), lambda);
      residual += c_abs2(c_sub(Hv_i, lv_i));
    }

    residual = sqrt(residual);

    char label[64];
    snprintf(label, sizeof label,
             "residual ||H v - lambda v|| for eigenpair %d", k);
    check_close(residual, 0.0, tol, label);
  }

  // Orthonormality check: V^dagger V ~= I
  for (int a = 0; a < n; a++) {
    for (int b = 0; b < n; b++) {
      complex_t dot = c_zero();
      for (int i = 0; i < n; i++) {
        dot = c_add(dot, c_mul(c_conj(CMAT(eig->eigenvectors, i, a)),
                               CMAT(eig->eigenvectors, i, b)));
      }
      double expected = (a == b) ? 1.0 : 0.0;
      double err = sqrt(c_abs2(c_sub(dot, c_real(expected))));

      char label[64];
      snprintf(label, sizeof label, "V^dagger V[%d][%d] identity", a, b);
      check_close(err, 0.0, tol, label);
    }
  }

  eigen_free(eig);
  cmatrix_free(H);
}

/* An 8-fold degenerate eigenvalue (+1) and a 2-fold one (-1) in a random
 * unitary basis */
static void test_large_degeneracy(void) {
  printf("  === Test: 8-fold degenerate 10x10 Hermitian (residual + "
         "orthonormality) ===\n");

  int n = 10;
  cmatrix_t *D = cmatrix_alloc(n, n);
  for (int i = 0; i < n * n; i++) {
    D->data[i] = c_zero();
  }

  const double eigvals_true[10] = {1, 1, 1, 1, 1, 1, 1, 1, -1, -1};
  for (int i = 0; i < n; i++) {
    CMAT(D, i, i) = c_real(eigvals_true[i]);
  }

  // NOTE: Random unitary-ish basis change via Gram-Schmidt on a random
  // complex matrix, to avoid testing only the trivial already-diagonal case
  cmatrix_t *Q = cmatrix_alloc(n, n);
  unsigned int seed = 42;

  for (int j = 0; j < n; j++) {
    complex_t col[10];

    for (int i = 0; i < n; i++) {
      seed = seed * 1103515245u + 12345u;
      double re = ((double)(seed % 10000) / 10000.0) - 0.5;

      seed = seed * 1103515245u + 12345u;
      double im = ((double)(seed % 10000) / 10000.0) - 0.5;

      col[i] = c_new(re, im);
    }

    for (int a = 0; a < j; a++) {
      complex_t dot = c_zero();

      for (int i = 0; i < n; i++) {
        dot = c_add(dot, c_mul(c_conj(CMAT(Q, i, a)), col[i]));
      }

      for (int i = 0; i < n; i++) {
        col[i] = c_sub(col[i], c_mul(dot, CMAT(Q, i, a)));
      }
    }

    double norm2 = 0.0;
    for (int i = 0; i < n; i++) {
      norm2 += c_abs2(col[i]);
    }

    double norm = sqrt(norm2);
    for (int i = 0; i < n; i++) {
      CMAT(Q, i, j) = c_scale(col[i], 1.0 / norm);
    }
  }

  // H = Q D Q^\dagger
  cmatrix_t *H = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      complex_t sum = c_zero();

      for (int k = 0; k < n; k++) {
        sum = c_add(sum, c_mul(c_scale(CMAT(Q, i, k), eigvals_true[k]),
                               c_conj(CMAT(Q, j, k))));
      }

      CMAT(H, i, j) = sum;
    }
  }

  cmatrix_free(Q);
  cmatrix_free(D);

  eigen_t *eig = cmatrix_eigh_complex(H);
  if (!eig) {
    check(0, "cmatrix_eigh_complex returned NULL");
    cmatrix_free(H);

    return;
  }

  double tol = 1e-8;
  double max_residual = 0.0;

  for (int k = 0; k < n; k++) {
    double lambda = eig->eigenvalues[k];
    double residual = 0.0;

    for (int i = 0; i < n; i++) {
      complex_t Hv_i = c_zero();

      for (int j = 0; j < n; j++) {
        Hv_i = c_add(Hv_i, c_mul(CMAT(H, i, j), CMAT(eig->eigenvectors, j, k)));
      }

      complex_t lv_i = c_scale(CMAT(eig->eigenvectors, i, k), lambda);
      residual += c_abs2(c_sub(Hv_i, lv_i));
    }

    residual = sqrt(residual);
    if (residual > max_residual) {
      max_residual = residual;
    }
  }

  printf("  max residual |H v - lambda v| = %.2e\n", max_residual);

  check_close(max_residual, 0.0, tol, "max eigenpair residual");

  double max_orth_err = 0.0;
  for (int a = 0; a < n; a++) {
    for (int b = 0; b < n; b++) {
      complex_t dot = c_zero();

      for (int i = 0; i < n; i++) {
        dot = c_add(dot, c_mul(c_conj(CMAT(eig->eigenvectors, i, a)),
                               CMAT(eig->eigenvectors, i, b)));
      }

      double expected = (a == b) ? 1.0 : 0.0;
      double err = sqrt(c_abs2(c_sub(dot, c_real(expected))));
      if (err > max_orth_err) {
        max_orth_err = err;
      }
    }
  }

  printf("  max orthonormality error |V^\\dagger V - I| = %.2e\n",
         max_orth_err);

  check_close(max_orth_err, 0.0, tol, "eigenbasis orthonormality");

  double eigvals_sorted[10];
  for (int i = 0; i < n; i++) {
    eigvals_sorted[i] = eig->eigenvalues[i];
  }

  for (int i = 1; i < n; i++) {
    double key = eigvals_sorted[i];
    int j = i - 1;
    while (j >= 0 && eigvals_sorted[j] > key) {
      eigvals_sorted[j + 1] = eigvals_sorted[j];
      j--;
    }

    eigvals_sorted[j + 1] = key;
  }

  const double eigvals_true_sorted[10] = {-1, -1, 1, 1, 1, 1, 1, 1, 1, 1};
  double max_eigval_err = 0.0;
  for (int i = 0; i < n; i++) {
    double err = fabs(eigvals_sorted[i] - eigvals_true_sorted[i]);
    if (err > max_eigval_err) {
      max_eigval_err = err;
    }
  }

  printf("  max eigenvalue error vs true {-1 x2, +1 x8} = %.2e\n",
         max_eigval_err);

  check_close(max_eigval_err, 0.0, tol, "eigenvalues match {-1 x2, +1 x8}");

  eigen_free(eig);
  cmatrix_free(H);
}

int main(void) {
  test_sigma_y();
  test_random_hermitian();
  test_large_degeneracy();

  if (failures == 0) {
    printf("\nAll test_complex_eigh checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_complex_eigh check(s) FAILED.\n", failures);
    return 1;
  }
}
