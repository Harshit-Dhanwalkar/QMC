# `anderson`

> Generated from `physics/anderson.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Anderson localisation on a 1D tight-binding chain with random on-site disorder.

**Category:** Condensed matter & many-body

```c
#include <qmc/physics/anderson.h>
```

**Start with:** [`anderson_disorder`](#anderson_disorder), [`anderson_evolve`](#anderson_evolve), [`anderson_norm`](#anderson_norm), [`anderson_width`](#anderson_width), [`anderson_ipr`](#anderson_ipr)

- Guide: [physics/anderson.md](../physics/anderson.md)
- Example: `examples/eg_79_anderson.c`
- Test (every call with expected results): `tests/test_anderson.c`
- Runtime help: `qmc_help(stdout, "anderson");` or `qmc-help anderson`

## Overview

```text
Anderson localisation on a 1D tight-binding chain (\hbar = 1)

  H = -t \sum_j (|j><j+1| + h.c.) + \sum_j \epsilon_j |j><j|

with n sites, open ends, hopping t > 0 and on-site energies \epsilon_j drawn
independently and uniformly from [-W/2, W/2] (W = disorder strength). Site j
sits at position x = j

NOTE: Without disorder (W = 0) a packet started on one site spreads
ballistically, \psi_j(\tau) = i^(j - j0) J_(j - j0)(2t\tau), so <x^2> = 2 t^2
\tau^2 exactly Any disorder, however weak, localises every eigenstate in 1D:
packet stops spreading and keeps an exponential profile \exp(-|x - x0|/xi)

For weak disorder inverse localisation length (Lyapunov exponent) of wave
amplitude at energy E inside band |E| < 2t is, to leading order,

    1/xi = W^2 / (24 (4 t^2 - E^2)),

which fails near band centre E = 0 (anomaly) and at band edges
```

## Functions (8)

[`anderson_disorder`](#anderson_disorder), [`anderson_evolve`](#anderson_evolve), [`anderson_norm`](#anderson_norm), [`anderson_width`](#anderson_width), [`anderson_ipr`](#anderson_ipr), [`anderson_energy`](#anderson_energy), [`anderson_lyapunov`](#anderson_lyapunov), [`anderson_lyapunov_weak`](#anderson_lyapunov_weak)

### `anderson_disorder`

```c
int anderson_disorder(double *eps, int n, double w, unsigned long long seed);
```

```text
Fill eps[0..n-1] with uniform samples in [-w/2, w/2] from a deterministic
generator seeded by `seed` (same seed -> same sequence on every platform)

Returns 0 on success, -1 on NULL eps, n < 1, or w < 0 / non-finite
```

### `anderson_evolve`

```c
int anderson_evolve(cvector_t *psi, const double *eps, double hop, double dt,
                    int steps);
```

```text
Advance psi (length n) by `steps` RK4 steps of size dt under chain
Hamiltonian above. eps has length n. Returns 0 on success, -1 on invalid
input (NULL pointers, n < 2, hop <= 0, dt <= 0, steps < 1, non-finite eps),
-2 on allocation failure. On failure psi is untouched

RK4 is not exactly unitary: with t = 1, dt = 0.01 norm drifts by about
5e-10 over 2000 steps for a clean chain and 2e-8 over 6000 steps at W = 6
```

### `anderson_norm`

```c
double anderson_norm(const cvector_t *psi);
```

```text
\sum_j |\psi_j|^2; NaN on NULL
```

### `anderson_width`

```c
double anderson_width(const cvector_t *psi);
```

```text
Root-mean-square width sqrt(<x^2> - <x>^2) in sites (x = j); NaN on NULL or
zero norm
```

### `anderson_ipr`

```c
double anderson_ipr(const cvector_t *psi);
```

```text
Inverse participation ratio sum |psi|^4 / (sum |psi|^2)^2. It is 1 for a
state on one site and 1/n for a state spread evenly over n sites; its
inverse counts how many sites state occupies. NaN on NULL or zero norm
```

### `anderson_energy`

```c
double anderson_energy(const cvector_t *psi, const double *eps, double hop);
```

```text
Energy expectation <psi|H|psi> / <psi|psi>; NaN on invalid input
```

### `anderson_lyapunov`

```c
double anderson_lyapunov(const double *eps, int n, double hop, double energy);
```

```text
Lyapunov exponent (inverse localisation length, per site) of amplitude at
energy E for given on-site energies, from transfer-matrix recursion
\psi_(j+1) = ((E - \eps_j)/t) \psi_j -  \psi_(j-1) with periodic
renormalisation: \gamma = (1/n) sum ln(growth). Self-averaging, so use a long
chain (>= 1e5 sites). NaN on invalid input
```

### `anderson_lyapunov_weak`

```c
double anderson_lyapunov_weak(double w, double hop, double energy);
```

```text
Weak-disorder estimate 1/xi = W^2 / (24 (4 t^2 - E^2)) for |E| < 2t; NaN
outside band or for invalid input
```
