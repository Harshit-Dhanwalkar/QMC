# `double_barrier`

> Generated from `physics/double_barrier.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Resonant tunnelling: transfer-matrix transmission through stacked square barriers, wavefunction profile and dwell time.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/double_barrier.h>
```

**Start with:** [`db_transmission`](#db_transmission), [`db_solve`](#db_solve), [`db_profile`](#db_profile), [`db_dwell_time`](#db_dwell_time)

- Example: `examples/eg_94_double_barrier.c`
- Test (every call with expected results): `tests/test_double_barrier.c`
- Runtime help: `qmc_help(stdout, "double_barrier");` or `qmc-help double_barrier`

## Overview

```text
Transmission through a stack of square barriers (transfer matrices), for
resonant tunnelling

Units \hbar^2 / 2m = 1, so E = k^2. A structure is n layers of potential v[i]
and thickness d[i] laid out from x = 0 to L = sum d; leads x < 0 and x > L
have V = 0. For a particle incident from left with unit amplitude transfer
matrix of (\psi, \psi') through a layer of local wave number q = sqrt(E - V)
is
[[\cos(qd), \sin(qd) / q], [-q \sin(qd), \cos(qd)]]
(\cosh/\sinh for E < V, [[1, d], [0, 1]] at E = V), a real unimodular matrix
and

    T = 4 / [ (M11 + M22)^2 + (k M12 - M21 / k)^2 ]

A single barrier gives T = 1 / (1 + V^2 \sinh^2(\kappa a) / (4 E (V - E)));
two identical barriers around a well show Fabry-Perot resonances with T = 1
at quasi-bound levels (n \pi / w)^2 of well, whose width shrinks as barriers
get opaque and whose dwell time grows as inverse width
```

## Types

### `?`

```c
enum { DB_MAX_LAYERS = 16 };
```

### `db_result_t`

```c
typedef struct {
  double transmission, reflection; /* |t|^2, |r|^2 */
  double r_re, r_im;               /* reflection amplitude at x = 0 */
} db_result_t;
```

## Functions (4)

[`db_transmission`](#db_transmission), [`db_solve`](#db_solve), [`db_profile`](#db_profile), [`db_dwell_time`](#db_dwell_time)

### `db_transmission`

```c
double db_transmission(int n_layers, const double *v, const double *d,
                       double energy);
```

```text
Transmission probability by closed transfer-matrix formula; NaN unless E
> 0, 1 <= n <= DB_MAX_LAYERS and d >= 0
```

### `db_solve`

```c
int db_solve(int n_layers, const double *v, const double *d, double energy,
             db_result_t *out);
```

```text
Complex reflection and transmission by matching plane waves; 0 on success, -1
on bad input
```

### `db_profile`

```c
int db_profile(int n_layers, const double *v, const double *d, double energy,
               double xmin, double xmax, int count, double *psi2);
```

```text
|\psi(x)|^2 for unit incident amplitude on `count` points from xmin to xmax
(leads included; structure spans 0 .. L). Returns 0, or -1 on bad input
(count >= 2, xmin < xmax)
```

### `db_dwell_time`

```c
double db_dwell_time(int n_layers, const double *v, const double *d,
                     double energy);
```

```text
Dwell time: integral of |\psi|^2 over structure divided by incident flux 2k;
NaN on bad input
```
