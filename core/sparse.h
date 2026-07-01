#ifndef QMC_SPARSE_H
#define QMC_SPARSE_H

#include "matrix.h"
#include "vector.h"

/* CSR (Compressed Sparse Row) format */
typedef struct {
  int nrows, ncols;
  int nnz;
  int *row_ptr;      /* size nrows+1 */
  int *col_ind;      /* size nnz */
  complex_t *values; /* size nnz */
} sparse_matrix_t;

/* Create/destroy */
sparse_matrix_t *sparse_alloc(int nrows, int ncols, int nnz);
void sparse_free(sparse_matrix_t *A);

/* Build from dense matrix (extract entries with |value| > tol) */
sparse_matrix_t *sparse_from_dense(const cmatrix_t *A, double tol);

/* Sparse matrix-vector multiply: y = A * x */
void sparse_mv(const sparse_matrix_t *A, const cvector_t *x, cvector_t *y);

/* For Hermitian (real symmetric) */
static inline void sparse_mv_hermitian(const sparse_matrix_t *A,
                                       const cvector_t *x, cvector_t *y) {
  sparse_mv(A, x, y);
}

/* Lanczos for lowest eigenvalues */
// HACK: stub for now */
typedef struct {
  int n;
  double *values;
  cmatrix_t *vectors;
} lanczos_result_t;

lanczos_result_t *lanczos_eigs(sparse_matrix_t *A, int k, int max_iter,
                               double tol);
void lanczos_free(lanczos_result_t *res);

#endif
