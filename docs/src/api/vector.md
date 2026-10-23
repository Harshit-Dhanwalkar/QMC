# `vector`

> Generated from `core/vector.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Complex vector type cvector_t: allocation, copy, dot product, norm and related operations.

**Category:** Core: data types & utilities

```c
#include <qmc/core/vector.h>
```

**Start with:** [`cvector_alloc`](#cvector_alloc), [`cvector_copy`](#cvector_copy), [`cvector_dot`](#cvector_dot), [`cvector_norm`](#cvector_norm), [`cvector_normalize`](#cvector_normalize)

- Runtime help: `qmc_help(stdout, "vector");` or `qmc-help vector`

## Types

### `cvector_t`

```c
typedef struct {
  complex_t *data;
  int n;
} cvector_t;
```

## Functions (10)

[`cvector_alloc`](#cvector_alloc), [`cvector_free`](#cvector_free), [`cvector_copy`](#cvector_copy), [`cvector_dot`](#cvector_dot), [`cvector_norm`](#cvector_norm), [`cvector_normalize`](#cvector_normalize), [`cvector_expect`](#cvector_expect), [`cvector_print`](#cvector_print), [`cvector_fill`](#cvector_fill), [`cvector_get`](#cvector_get)

### `cvector_alloc`

```c
cvector_t *cvector_alloc(int n);
```

```text
Allocate/Free
```

### `cvector_free`

```c
void cvector_free(cvector_t *v);
```

_No description in header yet._

### `cvector_copy`

```c
cvector_t *cvector_copy(const cvector_t *vec);
```

_No description in header yet._

### `cvector_dot`

```c
complex_t cvector_dot(const cvector_t *vec_a, const cvector_t *vec_b);
```

```text
Operations
```

### `cvector_norm`

```c
double cvector_norm(const cvector_t *vec);
```

_No description in header yet._

### `cvector_normalize`

```c
void cvector_normalize(cvector_t *vec);
```

_No description in header yet._

### `cvector_expect`

```c
double cvector_expect(const cvector_t *psi, const cvector_t *op_psi);
```

_No description in header yet._

### `cvector_print`

```c
void cvector_print(const cvector_t *vec, const char *label);
```

```text
Utilities
```

### `cvector_fill`

```c
void cvector_fill(cvector_t *vec, complex_t val);
```

_No description in header yet._

### `cvector_get`

```c
static inline complex_t cvector_get(const cvector_t *v, int index);
```

_No description in header yet._
