# `sparse`

> Generated from `core/sparse.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

CSR sparse matrices: construction from dense, matrix-vector products (including Hermitian).

**Category:** Core: data types & utilities

```c
#include <qmc/core/sparse.h>
```

**Start with:** [`sparse_alloc`](#sparse_alloc), [`sparse_from_dense`](#sparse_from_dense), [`sparse_mv`](#sparse_mv), [`sparse_mv_hermitian`](#sparse_mv_hermitian), [`lanczos_eigs`](#lanczos_eigs)

- Guide: [internals/sparse_solvers.md](../internals/sparse_solvers.md)
- Runtime help: `qmc_help(stdout, "sparse");` or `qmc-help sparse`

## Types

### `sparse_matrix_t`

```text
CSR (Compressed Sparse Row) format
```

```c
typedef struct {
  int nrows, ncols;
  int nnz;
  int *row_ptr;      /* size nrows+1 */
  int *col_ind;      /* size nnz */
  complex_t *values; /* size nnz */
} sparse_matrix_t;
```

### `lanczos_result_t`

```text
Lanczos for lowest eigenvalues
```

```c
typedef struct {
  int n; /* number of eigenpairs actually returned (== k passed to lanczos_eigs)
          */
  double *values;
  cmatrix_t *vectors;
} lanczos_result_t;
```

### `lanczos_tridiag_t`

```text
Lanczos tridiagonalization from a caller-supplied starting vector

Runs same three-term recurrence with full reorthogonalization starting from
v0, and returns raw tridiagonal coefficients (\alpha, \beta) (instead of
eigenpairs)

NOTE: v0 must already be normalized (||v0|| = 1); this is caller's
responsibility since normalization constant I0 = <v0|v0> (before normalizing)
is itself physically meaningful (e.g. total spectral weight / static
structure factor) and callers need it separately.

Returns NULL on invalid input or allocation failure. Otherwise returns a
lanczos_tridiag_t with m \alpha coefficients and m-1 \beta coefficients (m <=
max_iter, fewer if Krylov subspace collapses early)
```

```c
typedef struct {
  int m;
  double *alpha; /* size m */
  double *beta;  /* size m-1 (NULL if m == 1) */
} lanczos_tridiag_t;
```

## Functions (10)

[`sparse_alloc`](#sparse_alloc), [`sparse_free`](#sparse_free), [`sparse_from_dense`](#sparse_from_dense), [`sparse_mv`](#sparse_mv), [`sparse_mv_hermitian`](#sparse_mv_hermitian), [`lanczos_eigs`](#lanczos_eigs), [`lanczos_free`](#lanczos_free), [`lanczos_tridiagonalize`](#lanczos_tridiagonalize), [`lanczos_tridiag_free`](#lanczos_tridiag_free), [`cg_solve`](#cg_solve)

### `sparse_alloc`

```c
sparse_matrix_t *sparse_alloc(int nrows, int ncols, int nnz);
```

```text
Create/destroy
```

### `sparse_free`

```c
void sparse_free(sparse_matrix_t *sp_mat);
```

_No description in the header yet._

### `sparse_from_dense`

```c
sparse_matrix_t *sparse_from_dense(const cmatrix_t *sp_mat, double tol);
```

```text
Build from dense matrix (extract entries with |value| > tol)
```

### `sparse_mv`

```c
void sparse_mv(const sparse_matrix_t *sp_mat, const cvector_t *in_vec,
               cvector_t *out_vec);
```

```text
Sparse matrix-vector multiply: y = sp_mat * x
```

### `sparse_mv_hermitian`

```c
static inline void sparse_mv_hermitian(const sparse_matrix_t *sp_mat,
                                       const cvector_t *in_vec,
                                       cvector_t *out_vec);
```

```text
For Hermitian (real symmetric)
```

### `lanczos_eigs`

```c
lanczos_result_t *lanczos_eigs(const sparse_matrix_t *sp_mat, int k,
                               int max_iter, double tol);
```

```text
Lanczos iteration for k algebraically lowest eigenvalues/vectors of Hermitian
sparse matrix `sp_mat`, via a real tridiagonal Krylov projection with full
reorthogonalization

k       : number of lowest eigenvalues wanted.(1 <= k <= sp_mat->nrows)
max_iter: Lanczos steps to run (>= k); internally capped at sp_mat->nrows,
          since Krylov subspace can't exceed problem dimension
tol     : breakdown threshold for residual norm \beta_j; if \beta_j falls
          below this, invariant subspace found so far is used as-is

Returns NULL on invalid input, allocation failure, or if Krylov subspace
collapses (invariant subspace found) before k directions have been generated
Otherwise returns a lanczos_result_t with k lowest eigenvalues (ascending)
and their eigenvectors as columns of an n x k cmatrix_t
```

### `lanczos_free`

```c
void lanczos_free(lanczos_result_t *res);
```

_No description in the header yet._

### `lanczos_tridiagonalize`

```c
lanczos_tridiag_t *lanczos_tridiagonalize(const sparse_matrix_t *sp_mat,
                                          const cvector_t *v0, int max_iter,
                                          double tol);
```

_No description in the header yet._

### `lanczos_tridiag_free`

```c
void lanczos_tridiag_free(lanczos_tridiag_t *tridiag);
```

_No description in the header yet._

### `cg_solve`

```c
int cg_solve(const sparse_matrix_t *A, const cvector_t *b, cvector_t *x,
             int max_iter, double tol);
```

```text
Conjugate Gradient for Hermitian positive-definite sparse systems: solves
A x = b

A       : Hermitian positive-definite sparse matrix (not checked - a
          non-Hermitian or indefinite A will not converge, or will
          converge to a spurious solution)
b       : right-hand side (size A->nrows)
x       : in/out. On input, initial guess (zero vector is a fine default);
          on output, solution. Must be preallocated with
          x->n == b->n == A->nrows
max_iter: cap on CG iterations. In exact arithmetic, CG on an n x n  SPD
          system converges in at most n steps, but finite-precision CG can
          still make progress past that bound on ill-conditioned systems
          (loss of A-orthogonality between search directions), so this is not
          internally capped at A->nrows - pass a generous budget (e.g.
          several times n) for anything but a well-conditioned system
tol     : convergence threshold on residual norm  ||b - A x||

Returns 0 on convergence (residual norm < tol within max_iter iterations),
-1 on invalid input, allocation failure, or a breakdown (p^H A p == 0  before
convergence - can only happen if A is not actually positive definite)
```
