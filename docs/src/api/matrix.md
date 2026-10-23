# `matrix`

> Generated from `core/matrix.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Dense complex matrix type cmatrix_t: allocation, copy, multiply, transpose and related operations.

**Category:** Core: data types & utilities

```c
#include <qmc/core/matrix.h>
```

**Start with:** [`cmatrix_alloc`](#cmatrix_alloc), [`cmatrix_copy`](#cmatrix_copy), [`cmatrix_multiply`](#cmatrix_multiply), [`cmatrix_transpose`](#cmatrix_transpose), [`cmatrix_adjoint`](#cmatrix_adjoint)

- Test (every call with expected results): `tests/test_matrix.c`
- Runtime help: `qmc_help(stdout, "matrix");` or `qmc-help matrix`

## Types

### `cmatrix_t`

```c
typedef struct {
  complex_t *data; /* Row-major storage */
  int nrows, ncols;
} cmatrix_t;
```

### `eigen_t`

```text
Eigenvalues (for Hermitian matrices)
```

```c
typedef struct {
  int n;
  double *eigenvalues;
  cmatrix_t *eigenvectors;
} eigen_t;
```

## Macros and constants

```c
#define CMAT(m, i, j) ((m)->data[(i) * (m)->ncols + (j)])
```

```text
Access
```

## Functions (18)

[`cmatrix_alloc`](#cmatrix_alloc), [`cmatrix_free`](#cmatrix_free), [`cmatrix_copy`](#cmatrix_copy), [`cmatrix_multiply`](#cmatrix_multiply), [`cmatrix_transpose`](#cmatrix_transpose), [`cmatrix_adjoint`](#cmatrix_adjoint), [`cmatrix_scale`](#cmatrix_scale), [`cmatrix_mv`](#cmatrix_mv), [`cmatrix_add`](#cmatrix_add), [`cmatrix_kron`](#cmatrix_kron), [`cmatrix_lu_decomp`](#cmatrix_lu_decomp), [`cmatrix_solve`](#cmatrix_solve), [`cmatrix_eigh`](#cmatrix_eigh), [`eigen_free`](#eigen_free), [`cmatrix_set`](#cmatrix_set), [`cmatrix_get_row`](#cmatrix_get_row), [`cmatrix_permanent`](#cmatrix_permanent), [`cmatrix_print`](#cmatrix_print)

### `cmatrix_alloc`

```c
cmatrix_t *cmatrix_alloc(int nrows, int ncols);
```

```text
Allocate/Free
```

### `cmatrix_free`

```c
void cmatrix_free(cmatrix_t *mat);
```

_No description in header yet._

### `cmatrix_copy`

```c
cmatrix_t *cmatrix_copy(const cmatrix_t *mat);
```

_No description in header yet._

### `cmatrix_multiply`

```c
cmatrix_t *cmatrix_multiply(const cmatrix_t *left, const cmatrix_t *right);
```

```text
Operations
```

### `cmatrix_transpose`

```c
cmatrix_t *cmatrix_transpose(const cmatrix_t *matrix);
```

_No description in header yet._

### `cmatrix_adjoint`

```c
cmatrix_t *cmatrix_adjoint(const cmatrix_t *matrix);
```

_No description in header yet._

### `cmatrix_scale`

```c
void cmatrix_scale(cmatrix_t *matrix, complex_t scalar);
```

```text
Conjugate transpose
```

### `cmatrix_mv`

```c
cvector_t *cmatrix_mv(const cmatrix_t *mat,
                      const cvector_t *vec);
```

_No description in header yet._

### `cmatrix_add`

```c
cmatrix_t *cmatrix_add(const cmatrix_t *mat_a,
                       const cmatrix_t *mat_b);
```

```text
y = A * x
```

### `cmatrix_kron`

```c
cmatrix_t *cmatrix_kron(const cmatrix_t *mat_a,
                        const cmatrix_t *mat_b);
```

```text
c = a + b
```

### `cmatrix_lu_decomp`

```c
void cmatrix_lu_decomp(cmatrix_t *mat, int *pivot);
```

```text
Linear algebra
```

### `cmatrix_solve`

```c
cvector_t *cmatrix_solve(cmatrix_t *mat, const cvector_t *vec_b);
```

_No description in header yet._

### `cmatrix_eigh`

```c
eigen_t *cmatrix_eigh(const cmatrix_t *mat);
```

_No description in header yet._

### `eigen_free`

```c
void eigen_free(eigen_t *eigen_val);
```

```text
Hermitian eigendecomposition
```

### `cmatrix_set`

```c
static inline void cmatrix_set(cmatrix_t *m, int row, int col, complex_t val);
```

```text
Setter helper
```

### `cmatrix_get_row`

```c
cvector_t *cmatrix_get_row(const cmatrix_t *matrix, int row);
```

_No description in header yet._

### `cmatrix_permanent`

```c
complex_t cmatrix_permanent(const cmatrix_t *matrix);
```

_No description in header yet._

### `cmatrix_print`

```c
void cmatrix_print(const cmatrix_t *mat, const char *label);
```

```text
Utilities
```
