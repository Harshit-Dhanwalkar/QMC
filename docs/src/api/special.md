# `special`

> Generated from `core/special/special.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Special functions: Hermite, Laguerre, Legendre, Bessel, spherical harmonics and related polynomials.

**Category:** Core: special functions

```c
#include <qmc/core/special/special.h>
```

**Start with:** [`hermite`](#hermite), [`laguerre`](#laguerre), [`legendre`](#legendre)

- Guide: [internals/special.md](../internals/special.md)
- Test (every call with expected results): `tests/test_special.c`
- Runtime help: `qmc_help(stdout, "special");` or `qmc-help special`

## Functions (23)

[`hermite`](#hermite), [`hermite_array`](#hermite_array), [`hermite_deriv`](#hermite_deriv), [`hermite_zeros`](#hermite_zeros), [`hermite_zeros_all`](#hermite_zeros_all), [`hermite_function_stable`](#hermite_function_stable), [`laguerre`](#laguerre), [`laguerre_array`](#laguerre_array), [`legendre`](#legendre), [`assoc_legendre`](#assoc_legendre), [`legendre_array`](#legendre_array), [`spherical_harmonic_real`](#spherical_harmonic_real), [`spherical_harmonic_imag`](#spherical_harmonic_imag), [`spherical_harmonic`](#spherical_harmonic), [`sph_bessel_j`](#sph_bessel_j), [`sph_bessel_y`](#sph_bessel_y), [`sph_bessel_j_deriv`](#sph_bessel_j_deriv), [`sph_bessel_array`](#sph_bessel_array), [`riccati_bessel_j`](#riccati_bessel_j), [`riccati_bessel_y`](#riccati_bessel_y), [`riccati_bessel_j_deriv`](#riccati_bessel_j_deriv), [`factorial`](#factorial), [`log_factorial`](#log_factorial)

### `hermite`

```c
double hermite(int n, double x);
```

```text
Hermite polynomials
```

### `hermite_array`

```c
void hermite_array(int n, const double *x, int m, double *H);
```

_No description in header yet._

### `hermite_deriv`

```c
double hermite_deriv(int n, double x);
```

_No description in header yet._

### `hermite_zeros`

```c
double hermite_zeros(int n, int k);
```

_No description in header yet._

### `hermite_zeros_all`

```c
int hermite_zeros_all(int n, double *zeros);
```

_No description in header yet._

### `hermite_function_stable`

```c
double hermite_function_stable(int n, double xi);
```

```text
Numerically stable normalized Hermite function, via direct recurrence on
normalized function
```

### `laguerre`

```c
double laguerre(int n, double alpha, double x);
```

```text
Laguerre polynomials
```

### `laguerre_array`

```c
void laguerre_array(int n, double alpha, const double *x, int N, double *L);
```

_No description in header yet._

### `legendre`

```c
double legendre(int l, double x);
```

```text
Legendre polynomials
```

### `assoc_legendre`

```c
double assoc_legendre(int l, int m, double x);
```

_No description in header yet._

### `legendre_array`

```c
void legendre_array(int l, const double *x, int N, double *P);
```

_No description in header yet._

### `spherical_harmonic_real`

```c
double spherical_harmonic_real(int l, int m, double theta, double phi);
```

```text
Spherical harmonics
```

### `spherical_harmonic_imag`

```c
double spherical_harmonic_imag(int l, int m, double theta, double phi);
```

_No description in header yet._

### `spherical_harmonic`

```c
complex_t spherical_harmonic(int l, int m, double theta, double phi);
```

_No description in header yet._

### `sph_bessel_j`

```c
double sph_bessel_j(int l, double x);
```

```text
Spherical Bessel functions
```

### `sph_bessel_y`

```c
double sph_bessel_y(int l, double x);
```

_No description in header yet._

### `sph_bessel_j_deriv`

```c
double sph_bessel_j_deriv(int l, double x);
```

_No description in header yet._

### `sph_bessel_array`

```c
void sph_bessel_array(int lmax, double x, double *j, double *y);
```

_No description in header yet._

### `riccati_bessel_j`

```c
double riccati_bessel_j(int l, double x);
```

```text
Riccati-Bessel functions
```

### `riccati_bessel_y`

```c
double riccati_bessel_y(int l, double x);
```

_No description in header yet._

### `riccati_bessel_j_deriv`

```c
double riccati_bessel_j_deriv(int l, double x);
```

_No description in header yet._

### `factorial`

```c
double factorial(int n);
```

```text
Utilities
```

### `log_factorial`

```c
double log_factorial(int n);
```

_No description in header yet._
