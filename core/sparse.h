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
void sparse_free(sparse_matrix_t *sp_mat);

/* Build from dense matrix (extract entries with |value| > tol) */
sparse_matrix_t *sparse_from_dense(const cmatrix_t *sp_mat, double tol);

/* Sparse matrix-vector multiply: y = sp_mat * x */
void sparse_mv(const sparse_matrix_t *sp_mat, const cvector_t *in_vec,
               cvector_t *out_vec);

/* For Hermitian (real symmetric) */
static inline void sparse_mv_hermitian(const sparse_matrix_t *sp_mat,
                                       const cvector_t *in_vec,
                                       cvector_t *out_vec) {
  sparse_mv(sp_mat, in_vec, out_vec);
}

/* Lanczos for lowest eigenvalues */
typedef struct {
  int n;
  double *values;
  cmatrix_t *vectors;
} lanczos_result_t;

/*
 * Lanczos iteration for k algebraically lowest eigenvalues/vectors of Hermitian
 * sparse matrix `sp_mat`, via a real tridiagonal Krylov projection with full
 * reorthogonalization.
 *
 * k       : number of lowest eigenvalues wanted.(1 <= k <= sp_mat->nrows)
 * max_iter: Lanczos steps to run (>= k); internally capped at sp_mat->nrows, since
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
lanczos_result_t *lanczos_eigs(const sparse_matrix_t *sp_mat, int k,
                               int max_iter, double tol);
void lanczos_free(lanczos_result_t *res);

/* Lanczos tridiagonalization from a caller-supplied starting vector.
 *
 * Runs same three-term recurrence with full reorthogonalization starting from
 * v0, and returns raw tridiagonal coefficients (\alpha, \beta) (instead of
 * eigenpairs).
 *
 * NOTE: v0 must already be normalized (||v0|| = 1); this is caller's
 * responsibility since normalization constant I0 = <v0|v0> (before normalizing)
 * is itself physically meaningful (e.g. total spectral weight / static
 * structure factor) and callers need it separately.
 *
 * Returns NULL on invalid input or allocation failure. Otherwise returns a
 * lanczos_tridiag_t with m \alpha coefficients and m-1 \beta coefficients (m <=
 * max_iter, fewer if the Krylov subspace collapses early). Free with
 * lanczos_tridiag_free.
 */
typedef struct {
  int m;
  double *alpha; /* size m */
  double *beta;  /* size m-1 (NULL if m == 1) */
} lanczos_tridiag_t;

lanczos_tridiag_t *lanczos_tridiagonalize(const sparse_matrix_t *sp_mat,
                                          const cvector_t *v0, int max_iter,
                                          double tol);
void lanczos_tridiag_free(lanczos_tridiag_t *tridiag);

#endif
