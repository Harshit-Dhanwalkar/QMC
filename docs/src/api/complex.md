# `complex`

> Generated from `core/complex.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Complex number type complex_t and inline arithmetic (c_add, c_mul, c_exp, c_abs, ...).

**Category:** Core: data types & utilities

```c
#include <qmc/core/complex.h>
```

**Start with:** [`c_new`](#c_new), [`c_one`](#c_one), [`c_zero`](#c_zero), [`c_neg`](#c_neg), [`c_real`](#c_real)

- Test (every call with expected results): `tests/test_complex.c`
- Runtime help: `qmc_help(stdout, "complex");` or `qmc-help complex`

## Types

### `complex_t`

```c
typedef struct {
  double re, im;
} complex_t;
```

## Functions (18)

[`c_new`](#c_new), [`c_one`](#c_one), [`c_zero`](#c_zero), [`c_neg`](#c_neg), [`c_real`](#c_real), [`c_imag`](#c_imag), [`c_conj`](#c_conj), [`c_abs`](#c_abs), [`c_abs2`](#c_abs2), [`c_arg`](#c_arg), [`c_add`](#c_add), [`c_sub`](#c_sub), [`c_mul`](#c_mul), [`c_div`](#c_div), [`c_scale`](#c_scale), [`c_exp`](#c_exp), [`c_abs_sqr`](#c_abs_sqr), [`c_from_polar`](#c_from_polar)

### `c_new`

```c
static inline complex_t c_new(double re, double im);
```

```text
Basic operations
```

### `c_one`

```c
static inline complex_t c_one(void);
```

_No description in header yet._

### `c_zero`

```c
static inline complex_t c_zero(void);
```

_No description in header yet._

### `c_neg`

```c
static inline complex_t c_neg(complex_t z);
```

_No description in header yet._

### `c_real`

```c
static inline complex_t c_real(double x);
```

_No description in header yet._

### `c_imag`

```c
static inline complex_t c_imag(double y);
```

_No description in header yet._

### `c_conj`

```c
static inline complex_t c_conj(complex_t z);
```

_No description in header yet._

### `c_abs`

```c
static inline double c_abs(complex_t z);
```

```text
NOLINTNEXTLINE(clang-analyzer-core.UndefinedBinaryOperatorResult)
```

### `c_abs2`

```c
static inline double c_abs2(complex_t z);
```

_No description in header yet._

### `c_arg`

```c
static inline double c_arg(complex_t z);
```

_No description in header yet._

### `c_add`

```c
static inline complex_t c_add(complex_t a, complex_t b);
```

```text
Arithmetic
```

### `c_sub`

```c
static inline complex_t c_sub(complex_t a, complex_t b);
```

_No description in header yet._

### `c_mul`

```c
static inline complex_t c_mul(complex_t a, complex_t b);
```

_No description in header yet._

### `c_div`

```c
static inline complex_t c_div(complex_t a, complex_t b);
```

_No description in header yet._

### `c_scale`

```c
static inline complex_t c_scale(complex_t z, double s);
```

_No description in header yet._

### `c_exp`

```c
static inline complex_t c_exp(complex_t z);
```

_No description in header yet._

### `c_abs_sqr`

```c
static inline double c_abs_sqr(complex_t z);
```

```text
Polar coorindates
```

### `c_from_polar`

```c
static inline complex_t c_from_polar(double r, double theta);
```

_No description in header yet._
