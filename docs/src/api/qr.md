# `qr`

> Generated from `core/linalg/qr.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

QR decomposition with solve and inverse.

**Category:** Core: linear algebra

```c
#include <qmc/core/linalg/qr.h>
```

**Start with:** [`qr_decompose`](#qr_decompose), [`qr_solve`](#qr_solve), [`qr_invert`](#qr_invert)

- Guide: [internals/linalg.md](../internals/linalg.md)
- Runtime help: `qmc_help(stdout, "qr");` or `qmc-help qr`

## Functions (3)

[`qr_decompose`](#qr_decompose), [`qr_solve`](#qr_solve), [`qr_invert`](#qr_invert)

### `qr_decompose`

```c
int qr_decompose(cmatrix_t *A, cmatrix_t *Q);
```

```text
QR decomposition: A = Q R, where Q is unitary, R upper triangular.
A is overwritten with R (upper triangular part) and the Householder vectors
(lower part). Q is output as a separate matrix (must be preallocated or
NULL).

Returns 0 on success, -1 on failure.
```

### `qr_solve`

```c
int qr_solve(const cmatrix_t *A, const cmatrix_t *Q, const cvector_t *b,
             cvector_t *x);
```

```text
Solve least squares: min ||A x - b|| using QR.
A must have been QR decomposed (with qr_decompose).
b is RHS, x is solution (preallocated).

Returns 0 on success.
```

### `qr_invert`

```c
cmatrix_t *qr_invert(const cmatrix_t *A);
```

```text
Invert square matrix via QR (if non-singular).

Returns new matrix or NULL.
```
