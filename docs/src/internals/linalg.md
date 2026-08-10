# Linear Algebra Core

Dense and sparse linear algebra primitives underlying every physics module. Complex arithmetic (`core/complex.h`) and complex vectors/matrices (`core/vector.h`, `core/matrix.h`) are the base types; `core/linalg/` build decompositions and eigensolvers on top of them.

## Base types

```c
typedef struct { double re, im; } complex_t;

typedef struct {
  complex_t *data;
  int n;
} cvector_t;

typedef struct {
  complex_t *data; // Row-major storage
  int nrows, ncols;
} cmatrix_t;
#define CMAT(m, i, j) ((m)->data[(i) * (m)->ncols + (j)])

typedef struct {
  int n;
  double *eigenvalues;
  cmatrix_t *eigenvectors;
} eigen_t;
```

`complex_t` operations (`c_add`, `c_mul`, `c_div`, `c_conj`, `c_abs2`, `c_exp`, ...) are all `static inline` in `complex.h` - zero call overhead.

## Eigensolvers

QMC has three eigensolvers, chosen by matrix structure - using the wrong one is a real (and previously present) performance and correctness bug, not just a style choice:

| Function               | Matrix structure                    | Cost                                  | Header                   |
| ---------------------- | ----------------------------------- | ------------------------------------- | ------------------------ |
| `tridiag_eigh`         | real symmetric, tridiagonal         | $O(N^2)$                              | `linalg/tridiag_eigh.h`  |
| `cmatrix_eigh_complex` | general complex Hermitian           | via real-embedding                    | `linalg/complex_eigh.h`  |
| `cmatrix_eigh_generic` | general Hermitian (dense)           | LAPACK if available, else QR fallback | `linalg/eigen_generic.h` |
| `cmatrix_eigh`         | Hermitian (cyclic Jacobi rotations) | dense, $O(n^3)$ per full sweep        | `matrix.h`               |

```c
eigen_t *tridiag_eigh(const double *diag, const double *offdiag, int n);
```

Implicit QL algorithm with Wilkinson shift (EISPACK `tql2` / Numerical Recipes `tqli`). Takes the tridiagonal Hamiltonian as two flat arrays - no dense matrix is ever built. This is the correct solver for any finite-difference 1D (or radial) Schrödinger Hamiltonian, and replaced dense Jacobi diagonalization (`cmatrix_eigh_generic`) in several modules where the Hamiltonian is genuinely tridiagonal.

> **Known issue:** `tridiag_eigh` shows super-linear empirical scaling at larger grid sizes (suspected $O(N^{3-6})$ against the expected $O(N^2)$)

```c
eigen_t *cmatrix_eigh_complex(cmatrix_t *H);
```

For $H = A + iB$ where $H$ is genuinely complex-Hermitian (not just real tridiagonal), this builds the real symmetric $2N \times 2N$ embedding $\begin{pmatrix} A & -B \\ B & A \end{pmatrix}$ and diagonalizes that instead of implementing complex arithmetic directly in the eigensolver.

```c
eigen_t *cmatrix_eigh_generic(cmatrix_t *A);
```

General-purpose fallback: uses LAPACK if `USE_LAPACK` is defined at build time, otherwise falls back to `cmatrix_eigh` (QR-based).

## Direct solvers

```c
int *lu_decompose(cmatrix_t *A);                                // PA = LU
int lu_solve(const cmatrix_t *LU, const int *pivot,
             const cvector_t *b, cvector_t *x);
complex_t lu_det(const cmatrix_t *LU, const int *pivot);
cmatrix_t *lu_invert(const cmatrix_t *A);

int qr_decompose(cmatrix_t *A, cmatrix_t *Q);                   // A = QR
int qr_solve(const cmatrix_t *A, const cmatrix_t *Q,
             const cvector_t *b, cvector_t *x);
cmatrix_t *qr_invert(const cmatrix_t *A);

int svd_decompose(const cmatrix_t *A, cmatrix_t *U, cvector_t *S,
                  cmatrix_t *V);                                // A = U S V^\dagger
int svd_solve(const cmatrix_t *U, const cvector_t *S, const cmatrix_t *V,
             const cvector_t *b, cvector_t *x, double tol);
cmatrix_t *svd_pseudoinverse(const cmatrix_t *U, const cvector_t *S,
                             const cmatrix_t *V, double tol);
```

`lu_decompose` overwrites `A` in place with combined $L+U$ (unit diagonal not
stored). `qr_decompose` similarly overwrites `A` with $R$ and the Householder
vectors, with $Q$ returned separately.

## Sparse matrices

`core/sparse.h` provides CSR storage and matrix-vector products for large,
sparse Hamiltonians:

```c
typedef struct {
  int nrows, ncols, nnz;
  int *row_ptr;      // size nrows+1
  int *col_ind;      // size nnz
  complex_t *values; // size nnz
} sparse_matrix_t;

sparse_matrix_t *sparse_from_dense(const cmatrix_t *A, double tol);
void sparse_mv(const sparse_matrix_t *A, const cvector_t *x, cvector_t *y);
```

> **TODO:** `lanczos_eigs` (Lanczos iteration for the lowest few eigenvalues of a sparse matrix) is explicitly marked as a stub in the header (`// HACK: stub for now`) and not yet implemented.

## Utility access

```c
cvector_t *cvector_from_matrix_column(const cmatrix_t *m, int col);
```

Extracts a single eigenvector (or any matrix column) as a standalone `cvector_t` - the standard way to pull an individual eigenstate out of `eig->eigenvectors` after a call to any of the eigensolvers above.
