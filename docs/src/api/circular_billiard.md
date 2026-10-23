# `circular_billiard`

> Generated from `physics/circular_billiard.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Quantum billiard in a disk: Bessel functions and zeros, modes and Weyl counting.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/circular_billiard.h>
```

**Start with:** [`cb_bessel_j`](#cb_bessel_j), [`cb_bessel_zero`](#cb_bessel_zero), [`cb_modes`](#cb_modes), [`cb_mode_value`](#cb_mode_value)

- Example: `examples/eg_91_circular_billiard.c`
- Test (every call with expected results): `tests/test_circular_billiard.c`
- Runtime help: `qmc_help(stdout, "circular_billiard");` or `qmc-help circular_billiard`

## Overview

```text
Quantum particle in a circular well ("quantum corral" / circular drum)

Unit disc, hbar = 1, mass 1/2, so Schroedinger equation is
 -\nabla \psi = k^2 \psi with \psi = 0 on r = 1 and E = k^2

The eigenfunctions are:
     \psi(r, \theta) = N * J_m(k r) * \cos(m \theta)
  OR \psi(r, \theta) = N * J_m(k r) * \sin(m \theta)

with k = j_{m,n} n-th zero of Bessel function J_m

Every m >= 1 level is doubly degenerate (\cos and \sin); m = 0 is not
degenerate
A mode (m, n) has (n - 1) radial nodes and 2 m angular nodal lines
```

## Types

### `cb_mode_t`

```c
typedef struct {
  int m;       /* angular quantum number, 0, 1, 2, ... */
  int n;       /* radial quantum number, 1, 2, ... */
  int sine;    /* 0: \cos(m \theta), 1: \sin(m \theta) (only for m >= 1) */
  double k;    /* j_{m,n}; energy is k^2 */
  double norm; /* normalisation constant: integral of \psi^2 over disc is 1 */
} cb_mode_t;
```

## Functions (5)

[`cb_bessel_j`](#cb_bessel_j), [`cb_bessel_zero`](#cb_bessel_zero), [`cb_modes`](#cb_modes), [`cb_mode_value`](#cb_mode_value), [`cb_weyl_count`](#cb_weyl_count)

### `cb_bessel_j`

```c
double cb_bessel_j(int m, double x);
```

```text
Bessel function J_m(x) of integer order 0 <= m <= 200 for x >= 0; NaN
otherwise
```

### `cb_bessel_zero`

```c
double cb_bessel_zero(int m, int n);
```

```text
n-th positive zero of J_m (n >= 1), to about 1e-13; NaN on bad input
```

### `cb_modes`

```c
int cb_modes(double k_max, cb_mode_t *modes, int max_modes);
```

```text
All modes with k <= k_max, sorted by k (ties: cos before sin, smaller m
first) Writes at most max_modes entries and returns how many modes exist
(which may exceed max_modes)
Returns -1 on bad input
```

### `cb_mode_value`

```c
double cb_mode_value(const cb_mode_t *mode, double r, double theta);
```

```text
psi at radius r in [0, 1] and angle \theta
```

### `cb_weyl_count`

```c
double cb_weyl_count(double k);
```

```text
Weyl estimate of number of modes with k <= K: K^2 / 4 - K / 2
```
