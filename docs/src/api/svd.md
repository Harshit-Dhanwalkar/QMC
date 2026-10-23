# `svd`

> Generated from `core/linalg/svd.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Singular value decomposition, SVD solve and pseudoinverse.

**Category:** Core: linear algebra

```c
#include <qmc/core/linalg/svd.h>
```

**Start with:** [`svd_decompose`](#svd_decompose), [`svd_solve`](#svd_solve), [`svd_pseudoinverse`](#svd_pseudoinverse)

- Guide: [internals/linalg.md](../internals/linalg.md)
- Test (every call with expected results): `tests/test_svd.c`
- Runtime help: `qmc_help(stdout, "svd");` or `qmc-help svd`

## Functions (3)

[`svd_decompose`](#svd_decompose), [`svd_solve`](#svd_solve), [`svd_pseudoinverse`](#svd_pseudoinverse)

### `svd_decompose`

```c
int svd_decompose(const cmatrix_t *A, cmatrix_t *U, cvector_t *S, cmatrix_t *V);
```

```text
SVD: A = U * diag(S) * V^\dagger.
For complex A, we convert to real block matrix.
U and V are complex matrices, S is real vector.

Returns 0 on success.
```

### `svd_solve`

```c
int svd_solve(const cmatrix_t *U, const cvector_t *S, const cmatrix_t *V,
              const cvector_t *b, cvector_t *x, double tol);
```

```text
Solve least squares using SVD: min ||A x - b||.
SVD must have been computed.

Returns 0 on success.
```

### `svd_pseudoinverse`

```c
cmatrix_t *svd_pseudoinverse(const cmatrix_t *U, const cvector_t *S,
                             const cmatrix_t *V, double tol);
```

```text
Pseudoinverse: A+ = V * diag(1 / S) * U^\dagger.

Returns new matrix.
```
