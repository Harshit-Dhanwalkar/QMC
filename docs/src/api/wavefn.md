# `wavefn`

> Generated from `physics/wavefn.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Wavefunction container on a 1D grid: allocation, normalisation and probability density.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/wavefn.h>
```

**Start with:** [`wavefunction_alloc`](#wavefunction_alloc), [`wavefunction_normalize`](#wavefunction_normalize), [`wavefunction_free`](#wavefunction_free)

- Guide: [physics/wavefunctions.md](../physics/wavefunctions.md)
- Example: `examples/eg_05_uncertainty.c`
- Test (every call with expected results): `tests/test_wavefn.c`
- Runtime help: `qmc_help(stdout, "wavefn");` or `qmc-help wavefn`

## Types

### `wavefunction_t`

```c
typedef struct {
  cvector_t *psi; /* wavefunction values on grid */
  double *x;      /* position grid (owned) */ // NOTE: what is 'owned'?
  double dx;      /* grid spacing */
  int n;          /* number of grid points */
} wavefunction_t;
```

## Functions (15)

[`wavefunction_alloc`](#wavefunction_alloc), [`wavefunction_free`](#wavefunction_free), [`wavefunction_copy`](#wavefunction_copy), [`wavefunction_normalize`](#wavefunction_normalize), [`wavefunction_prob_density`](#wavefunction_prob_density), [`wavefunction_prob_in_interval`](#wavefunction_prob_in_interval), [`wavefunction_expect_x`](#wavefunction_expect_x), [`wavefunction_expect_x2`](#wavefunction_expect_x2), [`wavefunction_expect_p`](#wavefunction_expect_p), [`wavefunction_expect_p2`](#wavefunction_expect_p2), [`wavefunction_delta_x`](#wavefunction_delta_x), [`wavefunction_delta_p`](#wavefunction_delta_p), [`wavefunction_uncertainty_product`](#wavefunction_uncertainty_product), [`wavefunction_save`](#wavefunction_save), [`wavefunction_save_prob`](#wavefunction_save_prob)

### `wavefunction_alloc`

```c
wavefunction_t *wavefunction_alloc(int n);
```

```text
Allocate/free
```

### `wavefunction_free`

```c
void wavefunction_free(wavefunction_t *wavefn);
```

_No description in header yet._

### `wavefunction_copy`

```c
wavefunction_t *wavefunction_copy(const wavefunction_t *wavefn);
```

_No description in header yet._

### `wavefunction_normalize`

```c
void wavefunction_normalize(wavefunction_t *wavefn);
```

```text
Normalize (\int|\phi|^2 dx = 1)
```

### `wavefunction_prob_density`

```c
double *wavefunction_prob_density(const wavefunction_t *wavefn);
```

```text
Probability density array (|\phi|^2)
```

### `wavefunction_prob_in_interval`

```c
double wavefunction_prob_in_interval(const wavefunction_t *wavefn, double a,
                                     double b);
```

_No description in header yet._

### `wavefunction_expect_x`

```c
double wavefunction_expect_x(const wavefunction_t *wavefn);
```

```text
Expectation values (using central finite differences for p)
```

### `wavefunction_expect_x2`

```c
double wavefunction_expect_x2(const wavefunction_t *wavefn);
```

_No description in header yet._

### `wavefunction_expect_p`

```c
double wavefunction_expect_p(const wavefunction_t *wavefn);
```

_No description in header yet._

### `wavefunction_expect_p2`

```c
double wavefunction_expect_p2(const wavefunction_t *wavefn);
```

```text
uses FFT
```

### `wavefunction_delta_x`

```c
double wavefunction_delta_x(const wavefunction_t *wavefn);
```

```text
Uncertainty
```

### `wavefunction_delta_p`

```c
double wavefunction_delta_p(const wavefunction_t *wavefn);
```

_No description in header yet._

### `wavefunction_uncertainty_product`

```c
double wavefunction_uncertainty_product(const wavefunction_t *wavefn);
```

_No description in header yet._

### `wavefunction_save`

```c
void wavefunction_save(const wavefunction_t *wavefn, const char *filename);
```

```text
I/O
```

### `wavefunction_save_prob`

```c
void wavefunction_save_prob(const wavefunction_t *wavefn, const char *filename);
```

_No description in header yet._
