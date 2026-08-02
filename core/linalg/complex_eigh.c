/*
Complex-Hermitian eigensolver via real-embedding
*/

#include "complex_eigh.h"
#include "../complex.h"
#include "../matrix.h"
#include "eigen_generic.h"
#include <stdlib.h>

eigen_t *cmatrix_eigh_complex(cmatrix_t *H) {
  if (!H || H->nrows != H->ncols) {
    return NULL;
  }

  int n = H->nrows;
  int m2 = 2 * n;

  cmatrix_t *M = cmatrix_alloc(m2, m2);
  if (!M) {
    return NULL;
  }

  for (int i = 0; i < m2; i++) {
    for (int j = 0; j < m2; j++) {
      CMAT(M, i, j) = c_zero();
    }
  }

  // M = [[A,-B],[B,A]], H = A + iB
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      complex_t h = CMAT(H, i, j);

      CMAT(M, i, j) = c_real(h.re);
      CMAT(M, i, j + n) = c_real(-h.im);
      CMAT(M, i + n, j) = c_real(h.im);
      CMAT(M, i + n, j + n) = c_real(h.re);
    }
  }

  eigen_t *eig2n = cmatrix_eigh_generic(M);
  cmatrix_free(M);
  if (!eig2n) {
    return NULL;
  }
  if (eig2n->n != m2) {
    eigen_free(eig2n);
    return NULL;
  }

  eigen_t *result = malloc(sizeof *result);
  if (!result) {
    eigen_free(eig2n);
    return NULL;
  }

  result->n = n;
  result->eigenvalues = malloc(n * sizeof(double));
  result->eigenvectors = cmatrix_alloc(n, n);
  if (!result->eigenvalues || !result->eigenvectors) {
    free(result->eigenvalues);
    if (result->eigenvectors) {
      cmatrix_free(result->eigenvectors);
    }

    free(result);
    eigen_free(eig2n);

    return NULL;
  }

  for (int k = 0; k < n; k++) {
    int idx = 2 * k;
    double lambda = eig2n->eigenvalues[idx];
    if (idx + 1 < m2) {
      lambda = 0.5 * (lambda + eig2n->eigenvalues[idx + 1]);
    }
    result->eigenvalues[k] = lambda;

    for (int i = 0; i < n; i++) {
      double x = CMAT(eig2n->eigenvectors, i, idx).re;
      double y = CMAT(eig2n->eigenvectors, i + n, idx).re;
      CMAT(result->eigenvectors, i, k) = c_add(c_real(x), c_imag(y));
    }
  }

  eigen_free(eig2n);
  return result;
}
