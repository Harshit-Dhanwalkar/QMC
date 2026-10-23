# `eigen_generic`

> Generated from `core/linalg/eigen_generic.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Generic Hermitian eigensolver: LAPACK when built with USE_LAPACK, the built-in solver otherwise.

**Category:** Core: linear algebra

```c
#include <qmc/core/linalg/eigen_generic.h>
```

**Start with:** [`cmatrix_eigh_generic`](#cmatrix_eigh_generic)

- Guide: [internals/linalg.md](../internals/linalg.md)
- Runtime help: `qmc_help(stdout, "eigen_generic");` or `qmc-help eigen_generic`

## Functions (1)

[`cmatrix_eigh_generic`](#cmatrix_eigh_generic)

### `cmatrix_eigh_generic`

```c
eigen_t *cmatrix_eigh_generic(const cmatrix_t *A);
```

```text
Generic Hermitian eigen-solver.
Uses LAPACK if available (USE_LAPACK), otherwise falls back to
cmatrix_eigh (QR).
```
