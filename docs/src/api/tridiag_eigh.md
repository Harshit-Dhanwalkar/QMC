# `tridiag_eigh`

> Generated from `core/linalg/tridiag_eigh.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Eigenvalues and eigenvectors of real symmetric tridiagonal matrices (implicit QL).

**Category:** Core: linear algebra

```c
#include <qmc/core/linalg/tridiag_eigh.h>
```

**Start with:** [`tridiag_eigh`](#tridiag_eigh), [`tridiag_eigvals`](#tridiag_eigvals)

- Guide: [internals/linalg.md](../internals/linalg.md)
- Runtime help: `qmc_help(stdout, "tridiag_eigh");` or `qmc-help tridiag_eigh`

## Functions (2)

[`tridiag_eigh`](#tridiag_eigh), [`tridiag_eigvals`](#tridiag_eigvals)

### `tridiag_eigh`

```c
eigen_t *tridiag_eigh(const double *diag, const double *offdiag, int n);
```

```text
Eigen-decomposition of a real symmetric TRIDIAGONAL matrix via
implicit QL algorithm with Wilkinson shift (EISPACK tql2 /
Numerical-Recipes "tqli")

diag:     diagonal entries, size n
offdiag:  off-diagonal entries, size n-1;
          offdiag[i] connects row i and row i+1

Returns eigen_t* with eigenvalues ascending and eigenvectors as  columns of
an nxn cmatrix_t
```

### `tridiag_eigvals`

```c
eigen_t *tridiag_eigvals(const double *diag, const double *offdiag, int n);
```

```text
Eigenvalues only - same as tridiag_eigh, but skips all eigenvector
bookkeeping entirely (no working-array allocation, no per-rotation update, no
eigenvector permutation during sorting)

Returns eigen_t* with eigenvalues ascending and eigenvectors == NULL
Used instead of tridiag_eigh whenever caller only reads eig->eigenvalues
```
