# `lu`

> Generated from `core/linalg/lu.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

LU decomposition with partial pivoting: solve, determinant and inverse.

**Category:** Core: linear algebra

```c
#include <qmc/core/linalg/lu.h>
```

**Start with:** [`lu_decompose`](#lu_decompose), [`lu_solve`](#lu_solve), [`lu_det`](#lu_det), [`lu_invert`](#lu_invert)

- Guide: [internals/linalg.md](../internals/linalg.md)
- Runtime help: `qmc_help(stdout, "lu");` or `qmc-help lu`

## Functions (4)

[`lu_decompose`](#lu_decompose), [`lu_solve`](#lu_solve), [`lu_det`](#lu_det), [`lu_invert`](#lu_invert)

### `lu_decompose`

```c
int *lu_decompose(cmatrix_t *A);
```

```text
LU decomposition with partial pivoting: PA = LU.

 Returns pivot array (int*) of size n, or NULL on failure.
 The input matrix A is overwritten with L+U (unit diagonal not stored).
 The pivot array indicates row permutations.
```

### `lu_solve`

```c
int lu_solve(const cmatrix_t *LU, const int *pivot, const cvector_t *b,
             cvector_t *x);
```

```text
Solve A x = b using LU decomposition.
LU and pivot must come from lu_decompose.
b is the right-hand side; x is the solution (must be preallocated).

Returns 0 on success, -1 on failure.
```

### `lu_det`

```c
complex_t lu_det(const cmatrix_t *LU, const int *pivot);
```

```text
Determinant of A given LU and pivot.

Returns complex determinant.
```

### `lu_invert`

```c
cmatrix_t *lu_invert(const cmatrix_t *A);
```

```text
Inverse of A using LU decomposition.

 Returns new cmatrix_t* or NULL on failure.
```
