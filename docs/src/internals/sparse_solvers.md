# Sparse Iterative Solvers

`core/sparse.h`/`core/sparse.c` and `core/shift_invert.h`/`.c` provide three
iterative solvers built on top of the library's CSR sparse-matrix
representation, used wherever a Hilbert space is too large to diagonalize
densely - most directly in [DMRG](../physics/dmrg.md) and other
large-lattice methods.

## Sparse Matrix Storage

```c
typedef struct {
  int nrows, ncols;
  int nnz;
  int *row_ptr;      /* size nrows+1 */
  int *col_ind;      /* size nnz */
  complex_t *values; /* size nnz */
} sparse_matrix_t;

sparse_matrix_t *sparse_alloc(int nrows, int ncols, int nnz);
void sparse_free(sparse_matrix_t *sp_mat);
sparse_matrix_t *sparse_from_dense(const cmatrix_t *sp_mat, double tol);
void sparse_mv(const sparse_matrix_t *sp_mat, const cvector_t *in_vec,
               cvector_t *out_vec);
```

Standard Compressed Sparse Row (CSR) format. `sparse_from_dense` extracts
entries with `|value| > tol` from a dense `cmatrix_t` - the easiest way to
build a `sparse_matrix_t` for testing, or for small systems where dense
construction is convenient even though the solver itself is sparse.

All three solvers below only ever touch `A` through `sparse_mv` (`y = A x`),
never `A`'s entries directly, so none of them assume anything about how `A`
was built.

## Lanczos: Extremal Eigenvalues

```c
typedef struct {
  int n;              /* number of eigenpairs actually returned (== k) */
  double *values;      /* ascending */
  cmatrix_t *vectors;  /* n_A x k, eigenvectors as columns */
} lanczos_result_t;

lanczos_result_t *lanczos_eigs(const sparse_matrix_t *sp_mat, int k,
                               int max_iter, double tol);
void lanczos_free(lanczos_result_t *res);
```

Finds the `k` algebraically lowest eigenvalues/eigenvectors of a Hermitian
sparse matrix via Lanczos tridiagonalization (a three-term recurrence
building a Krylov subspace) with full reorthogonalization, followed by
diagonalizing the resulting real tridiagonal matrix. `max_iter` is
internally capped at `sp_mat->nrows`, since the Krylov subspace can't exceed
the problem dimension; `tol` is a breakdown threshold on the residual norm
$\beta_j$ - if it falls below `tol`, an invariant subspace has been found and
is used as-is rather than continuing.

Cheap per iteration (one `sparse_mv` plus $O(n)$ vector work), but **only
reaches the extremal end of the spectrum** - it has no way to target
eigenvalues in the interior. For those, see `shift_invert_eigs` below.

```c
sparse_matrix_t *H = /* Hermitian sparse Hamiltonian */;
lanczos_result_t *res = lanczos_eigs(H, 3, 200, 1e-10); // 3 lowest eigenvalues
if (res) {
  for (int i = 0; i < res->n; i++) {
    printf("E_%d = %.6f\n", i, res->values[i]);
  }
  lanczos_free(res);
}
```

There's also a lower-level entry point for when you need the raw
tridiagonal coefficients rather than eigenpairs - e.g. to compute a spectral
function or static structure factor from a physically meaningful starting
vector $v_0$ (rather than a random one):

```c
typedef struct {
  int m;
  double *alpha; /* size m */
  double *beta;  /* size m-1 (NULL if m == 1) */
} lanczos_tridiag_t;

lanczos_tridiag_t *lanczos_tridiagonalize(const sparse_matrix_t *sp_mat,
                                          const cvector_t *v0, int max_iter,
                                          double tol);
void lanczos_tridiag_free(lanczos_tridiag_t *tridiag);
```

`v0` must already be normalized ($\lVert v_0 \rVert = 1$) - normalizing it
internally would throw away $I_0 = \langle v_0 | v_0 \rangle$ (the
pre-normalization norm), which callers doing spectral-function-type
calculations need separately as the total spectral weight.

## Conjugate Gradient: Linear Systems

```c
int cg_solve(const sparse_matrix_t *A, const cvector_t *b, cvector_t *x,
             int max_iter, double tol);
```

Solves $Ax = b$ for Hermitian positive-definite `A` via the standard
Conjugate Gradient recurrence (not checked - a non-Hermitian or indefinite
`A` will not converge, or will converge to a spurious solution). `x` is
both the initial guess on input (the zero vector is a fine default) and the
solution on output; `tol` is a convergence threshold on the residual norm
$\lVert b - Ax \rVert$.

Unlike `lanczos_eigs`, `max_iter` is **not** capped at `A`'s dimension: in
exact arithmetic CG converges in at most $n$ steps, but finite-precision CG
on an ill-conditioned system can still make real progress past that
theoretical bound (loss of $A$-orthogonality between search directions
accumulates rounding error that more iterations can work back down), so the
caller's `max_iter` is honored directly - pass a generous budget (several
times $n$) for anything but a well-conditioned system.

Returns `0` on convergence, `-1` on invalid input, allocation failure, or a
breakdown ($p^H A p = 0$ before convergence, which can only happen if `A`
isn't actually positive definite).

```c
cvector_t *x = cvector_alloc(n);
cvector_fill(x, c_zero()); // zero initial guess
if (cg_solve(A, b, x, 500, 1e-10) == 0) {
  /* x holds the solution */
}
```

## Shift-Invert: Interior Eigenvalues

```c
typedef struct {
  int n;
  double *values;      /* eigenvalues nearest sigma, ascending */
  cmatrix_t *vectors;  /* n_A x k, eigenvectors as columns */
} shift_invert_result_t;

shift_invert_result_t *shift_invert_eigs(const sparse_matrix_t *A,
                                         double sigma, int k, int max_iter,
                                         double tol);
void shift_invert_result_free(shift_invert_result_t *res);
```

Finds the `k` eigenvalues nearest an arbitrary shift `sigma`, via inverse
iteration on $(A - \sigma I)$: a single dense LU factorization of the
shifted matrix is computed once up front, then reused across both all `k`
eigenvectors and every iteration of each. This is what makes it possible to
target interior eigenvalues at all - repeatedly solving $(A - \sigma I) w =
v$ and renormalizing converges $w$ toward whichever eigenvector of $A$ is
closest to $\sigma$, far faster than $A$'s other eigenvectors, since
inversion amplifies small eigenvalues of $(A - \sigma I)$ (i.e. eigenvalues
of $A$ near $\sigma$) the most.

For `k > 1`, each subsequent eigenvector is deflated (Gram-Schmidt
orthogonalized) against all previously-found ones at every iteration, so
the solver converges to progressively farther eigenvalues from `sigma`
rather than repeating the same one. `tol` is a convergence threshold on the
overlap between successive iterates, $|1 - |\langle v_{\text{prev}} | v
\rangle||$.

Because it needs one $O(n^3)$ dense factorization up front, this is meant
for matrices small/dense enough that a $n \times n$ dense factorization is
affordable - not a substitute for `lanczos_eigs` on very large sparse
systems, but the only option here for eigenvalues away from the spectrum's
edges.

Returns `NULL` on invalid input, allocation failure, or if $(A - \sigma I)$
is exactly singular - i.e. `sigma` is exactly an eigenvalue of `A`, where
shift-invert is undefined; perturb `sigma` slightly and retry.

```c
sparse_matrix_t *H = /* Hermitian sparse Hamiltonian, small enough to
                       * densify for one LU factorization */;

/* 2 eigenvalues nearest E=1.5, e.g. states around a known resonance */
shift_invert_result_t *res = shift_invert_eigs(H, 1.5, 2, 100, 1e-10);
if (res) {
  for (int i = 0; i < res->n; i++) {
    printf("E_%d = %.6f\n", i, res->values[i]);
  }

  shift_invert_result_free(res);
}
```
