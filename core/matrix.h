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
void cmatrix_free(cmatrix_t *mat);
cmatrix_t *cmatrix_copy(const cmatrix_t *mat);

/* Access */
#define CMAT(m, i, j) ((m)->data[(i) * (m)->ncols + (j)])

/* Operations */
cmatrix_t *cmatrix_multiply(const cmatrix_t *left, const cmatrix_t *right);
cmatrix_t *cmatrix_transpose(const cmatrix_t *matrix);
cmatrix_t *cmatrix_adjoint(const cmatrix_t *matrix); /* Conjugate transpose */
void cmatrix_scale(cmatrix_t *matrix, complex_t scalar);
cvector_t *cmatrix_mv(const cmatrix_t *mat,
                      const cvector_t *vec); /* y = A * x */
cmatrix_t *cmatrix_add(const cmatrix_t *mat_a,
                       const cmatrix_t *mat_b); /* c = a + b */

/* Linear algebra */
void cmatrix_lu_decomp(cmatrix_t *mat, int *pivot);
cvector_t *cmatrix_solve(cmatrix_t *mat, const cvector_t *vec_b);

/* Eigenvalues (for Hermitian matrices) */
typedef struct {
  int n;
  double *eigenvalues;
  cmatrix_t *eigenvectors;
} eigen_t;

eigen_t *cmatrix_eigh(const cmatrix_t *mat); /* Hermitian eigendecomposition */
void eigen_free(eigen_t *eigen_val);

/* Setter helper */
static inline void cmatrix_set(cmatrix_t *m, int row, int col, complex_t val) {
  m->data[row * m->ncols + col] = val;
}
cvector_t *cmatrix_get_row(const cmatrix_t *matrix, int row);
complex_t cmatrix_permanent(const cmatrix_t *matrix);

/* Utilities */
void cmatrix_print(const cmatrix_t *mat, const char *label);

#endif
