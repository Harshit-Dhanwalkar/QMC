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
typedef struct {
  int n;
  double *values;
  cmatrix_t *vectors;
} lanczos_result_t;

/*
 * Lanczos iteration for k algebraically lowest eigenvalues/vectors of Hermitian
 * sparse matrix A, via a real tridiagonal Krylov projection with full
 * reorthogonalization.
 *
 * k       : number of lowest eigenvalues wanted.(1 <= k <= A->nrows)
 * max_iter: Lanczos steps to run (>= k); internally capped at A->nrows, since
 *           the Krylov subspace can't exceed problem dimension.
 * tol     : breakdown threshold for residual norm \beta_j; if \beta_j falls
 *           below this, invariant subspace found so far is used as-is.
 *
 * Returns NULL on invalid input, allocation failure, or if Krylov subspace
 * collapses (invariant subspace found) before k directions have been generated.
 * Otherwise returns a lanczos_result_t with k lowest eigenvalues (ascending)
 * and their eigenvectors as columns of an n x k cmatrix_t. Free with
 * lanczos_free.
 */
lanczos_result_t *lanczos_eigs(const sparse_matrix_t *A, int k, int max_iter,
                               double tol);
void lanczos_free(lanczos_result_t *res);

#endif
