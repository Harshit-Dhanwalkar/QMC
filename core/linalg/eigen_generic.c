#include "eigen_generic.h"
#include "../matrix.h"
#include <stdio.h>
#include <stdlib.h>

// Forward declare the QR-based solver
extern eigen_t *cmatrix_eigh(const cmatrix_t *A);

#ifdef USE_LAPACK
#include <lapacke.h>

/*
 * Generic (real or complex) Hermitian solver via LAPACK's zheev, called on
 * every cmatrix_t regardless of whether its imaginary parts happen to be zero.
 */
/* NOTE: Unlike a real-only dsyev path, this can never silently drop a complex
 * Hermitian matrix's imaginary part. complex_t's layout ({double re, im;}) is
 * bit-compatible with LAPACK's lapack_complex_double (two consecutive doubles),
 * so cmatrix_t data can be reinterpreted in place with no copy/conversion
 * beyond what zheev itself needs (it overwrites its input array with the
 * eigenvectors, so a copy of the input is still required to preserve the
 * aller's A).
 */

eigen_t *cmatrix_eigh_lapack(const cmatrix_t *A) {
  if (!A || A->nrows != A->ncols) {
    return NULL;
  }

  int n = A->nrows;
  eigen_t *result = malloc(sizeof(eigen_t));
  if (!result) {
    return NULL;
  }

  result->n = n;
  result->eigenvalues = malloc((size_t)n * sizeof(double));
  result->eigenvectors = cmatrix_alloc(n, n);
  if (!result->eigenvalues || !result->eigenvectors) {
    free(result->eigenvalues);

    if (result->eigenvectors) {
      cmatrix_free(result->eigenvectors);
    }

    free(result);

    return NULL;
  }

  /* NOTE: zheev overwrites its input with the eigenvectors, so copy A into the
   * (row-major) buffer we're about to hand it : CMAT's row-major layout
   * matches LAPACK_ROW_MAJOR directly, no transpose needed.
   */
  for (int i = 0; i < n * n; i++) {
    result->eigenvectors->data[i] = A->data[i];
  }

  int info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', n,
                           (lapack_complex_double *)result->eigenvectors->data,
                           n, result->eigenvalues);
  if (info != 0) {
    fprintf(stderr, "LAPACK zheev failed with error %d\n", info);

    cmatrix_free(result->eigenvectors);
    free(result->eigenvalues);
    free(result);

    return NULL;
  }

  return result;
}
#else

// Fallback when LAPACK is not available - use the built‑in QR eigensolver.
eigen_t *cmatrix_eigh_lapack(const cmatrix_t *A) {
  return cmatrix_eigh(A);
}

#endif

eigen_t *cmatrix_eigh_generic(const cmatrix_t *A) {
  return cmatrix_eigh_lapack(A);
}
