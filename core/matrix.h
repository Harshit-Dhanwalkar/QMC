#ifndef QMC_MATRIX_H
#define QMC_MATRIX_H

#include "complex.h"
#include "vector.h"

typedef struct {
  complex_t *data; /* Row-major storage */
  int nrows, ncols;
} cmatrix_t;

/* Allocate/Free */
cmatrix_t *cmatrix_alloc(int nrows, int ncols);
void cmatrix_free(cmatrix_t *m);
cmatrix_t *cmatrix_copy(const cmatrix_t *m);

/* Access */
#define CMAT(m, i, j) ((m)->data[(i) * (m)->ncols + (j)])

/* Operations */
cmatrix_t *cmatrix_multiply(const cmatrix_t *a, const cmatrix_t *b);
cmatrix_t *cmatrix_transpose(const cmatrix_t *m);
cmatrix_t *cmatrix_adjoint(const cmatrix_t *m); /* Conjugate transpose */
void cmatrix_scale(cmatrix_t *m, complex_t s);

/* Linear algebra */
void cmatrix_lu_decomp(cmatrix_t *a, int *pivot);
cvector_t *cmatrix_solve(cmatrix_t *a, const cvector_t *b);

/* Eigenvalues (for Hermitian matrices) */
typedef struct {
  int n;
  double *eigenvalues;
  cmatrix_t *eigenvectors;
} eigen_t;

eigen_t *cmatrix_eigh(cmatrix_t *a); /* Hermitian eigendecomposition */
void eigen_free(eigen_t *e);

/* Utilities */
void cmatrix_print(const cmatrix_t *m, const char *label);

#endif
