// TODO:

// eigen_generic.c
#include "eigen_generic.h"
#include "../matrix.h"
#include <stdio.h>
#include <stdlib.h>

// Forward declare the QR-based solver
extern eigen_t *cmatrix_eigh(cmatrix_t *A);

// FIX: later, allocate matrix and copy result from LAPACK's column‑major
// storage
#ifdef USE_LAPACK
#include <lapacke.h>

eigen_t *cmatrix_eigh_lapack(cmatrix_t *A) {
  if (!A || A->nrows != A->ncols)
    return NULL;
  int n = A->nrows;
  double *mat = malloc(n * n * sizeof(double));
  double *eigvals = malloc(n * sizeof(double));
  if (!mat || !eigvals) {
    free(mat);
    free(eigvals);
    return NULL;
  }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      mat[j * n + i] = CMAT(A, i, j).re; // assumes real symmetric

  int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, mat, n, eigvals);
  if (info != 0) {
    fprintf(stderr, "LAPACK dsyev failed with error %d\n", info);
    free(mat);
    free(eigvals);
    return NULL;
  }

  eigen_t *result = malloc(sizeof(eigen_t));
  if (!result) {
    free(mat);
    free(eigvals);
    return NULL;
  }

  result->n = n;
  result->eigenvalues = eigvals;
  result->eigenvectors = malloc(n * sizeof(cvector_t));
  for (int i = 0; i < n; i++) {
    result->eigenvectors[i].n = n;
    result->eigenvectors[i].data = malloc(n * sizeof(complex_t));
    for (int j = 0; j < n; j++)
      result->eigenvectors[i].data[j] = c_real(mat[i * n + j]);
  }

  free(mat);

  return result;
}
#endif

eigen_t *cmatrix_eigh_generic(cmatrix_t *A) {
#ifdef USE_LAPACK
  return cmatrix_eigh_lapack(A);
#else
  return cmatrix_eigh(A); // existing QR implementation;
#endif
}
