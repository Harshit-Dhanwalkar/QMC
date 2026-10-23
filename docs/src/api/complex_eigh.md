# `complex_eigh`

> Generated from `core/linalg/complex_eigh.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Eigen-decomposition of complex Hermitian matrices via the real-embedding trick (cmatrix_eigh_complex).

**Category:** Core: linear algebra

```c
#include <qmc/core/linalg/complex_eigh.h>
```

**Start with:** [`cmatrix_eigh_complex`](#cmatrix_eigh_complex)

- Guide: [internals/linalg.md](../internals/linalg.md)
- Test (every call with expected results): `tests/test_complex_eigh.c`
- Runtime help: `qmc_help(stdout, "complex_eigh");` or `qmc-help complex_eigh`

## Functions (1)

[`cmatrix_eigh_complex`](#cmatrix_eigh_complex)

### `cmatrix_eigh_complex`

```c
eigen_t *cmatrix_eigh_complex(cmatrix_t *H);
```

```text
Eigen-decomposition of general complex Hermitian matrix via real-embedding
trick.

For H = A + iB, build real symmetric 2N x 2N matrix
  M = [[ A, -B ],
       [ B,  A ]]
```
