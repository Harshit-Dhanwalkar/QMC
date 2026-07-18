/*
Test: complex-Hermitian eigensolver (cmatrix_eigh_complex).

1. 2x2: Pauli \sigma_y matrix has known eigenvalues \pm 1 and known eigenvectors
(1,+-i)/\sqrt2
2. General NxN Hermitian check:
   verify ||H v_k - lambda_k v_k|| is small for every eigenpair, eigenvectors
  are orthonormal (V^\dagger V ~= I). This checks defining property of an
  eigendecomposition directly rather than needing a closed-form reference for
  arbitrary matrix.
*/

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/angular.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int test_sigma_y(void) {
  cmatrix_t *H = cmatrix_alloc(2, 2);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      CMAT(H, i, j) = sigma_y[i * 2 + j];

  eigen_t *eig = cmatrix_eigh_complex(H);
  cmatrix_free(H);
  if (!eig) {
    printf("  FAIL: cmatrix_eigh_complex returned NULL\n");
    return 1;
  }

  int fail = 0;
  double expected[2] = {-1.0, 1.0};
  for (int k = 0; k < 2; k++) {
    double err = fabs(eig->eigenvalues[k] - expected[k]);
    printf("  eigenvalue[%d]: got=%.6f expected=%.6f err=%.2e\n", k,
           eig->eigenvalues[k], expected[k], err);
    if (err > 1e-8)
      fail = 1;
  }

  // Eigenvector check up to phase: |v . conj(v_expected)| should be ~1.
  // v_expected for lambda=-1 is (1,-i)/\sqrt(2), for \lambda=+1 is
  // (1,i)/\sqrt(2).
  double inv_sqrt2 = 1.0 / sqrt(2.0);
  complex_t v_exp[2][2] = {{c_real(inv_sqrt2), c_imag(-inv_sqrt2)},
                           {c_real(inv_sqrt2), c_imag(inv_sqrt2)}};
  for (int k = 0; k < 2; k++) {
    complex_t overlap = c_zero();
    for (int i = 0; i < 2; i++) {
      complex_t vi = CMAT(eig->eigenvectors, i, k);
      overlap = c_add(overlap, c_mul(c_conj(vi), v_exp[k][i]));
    }

    double mag = sqrt(c_abs2(overlap));
    printf("  eigenvector[%d] overlap magnitude: %.6f (expected ~1)\n", k, mag);
    if (fabs(mag - 1.0) > 1e-6)
      fail = 1;
  }

  eigen_free(eig);
  return fail;
}

static int test_random_hermitian(void) {
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
    printf("  FAIL: cmatrix_eigh_complex returned NULL\n");
    cmatrix_free(H);
    return 1;
  }

  int fail = 0;
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
      complex_t diff = c_sub(Hv_i, lv_i);
      residual += c_abs2(diff);
    }

    residual = sqrt(residual);
    printf("  eigenpair %d: \\lambda=%.6f residual=%.2e\n", k, lambda,
           residual);

    if (residual > tol)
      fail = 1;
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
      if (err > tol) {
        printf("  FAIL: V^\\dagger V[%d][%d] err=%.2e\n", a, b, err);
        fail = 1;
      }
    }
  }

  eigen_free(eig);
  cmatrix_free(H);
  return fail;
}

int main(void) {
  int failed = 0;

  printf("Pauli \\sigma_y (known eigenvalues +-1):\n");
  failed += test_sigma_y();

  printf("Random 5x5 Hermitian (residual + orthonormality check):\n");
  failed += test_random_hermitian();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }

  printf("PASS\n");
  return 0;
}
