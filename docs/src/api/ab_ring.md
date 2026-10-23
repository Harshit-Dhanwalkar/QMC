# `ab_ring`

> Generated from `physics/ab_ring.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Aharonov-Bohm ring: energy levels versus flux, persistent current and flux-dependent transmission.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/ab_ring.h>
```

**Start with:** [`ar_level`](#ar_level), [`ar_levels_sorted`](#ar_levels_sorted), [`ar_energy`](#ar_energy), [`ar_current`](#ar_current), [`ar_scatter`](#ar_scatter)

- Example: `examples/eg_93_ab_ring.c`
- Test (every call with expected results): `tests/test_ab_ring.c`
- Runtime help: `qmc_help(stdout, "ab_ring");` or `qmc-help ab_ring`

## Overview

```text
Aharonov-Bohm ring: a tight-binding ring of N sites with hopping t threaded
by magnetic flux \Phi.

Gauge: flux phase exp(2 \pi i \phi), \phi = \Phi / \Phi_0, sits on one bond;
eigenstates are plane waves with energies E_n = -2 t \cos(2 pi (n - \phi) /
N), n = 0 .. N-1, periodic in phi with period 1

 - Persistent current of ne spinless electrons:
    I = -dE_total/d\phi, with E_total sum of ne lowest levels
   It is a periodic sawtooth-like function of flux, diamagnetic or
   paramagnetic depending on parity of ne
   For odd ne = 2m+1 and |\phi| < 1/2:
    I = -(4 \pi t / N) \sin(2 \pi \phi / N) \sin((2m+1) \pi / N) /
         \sin(\pi / N)
   for even ne = 2m and 0 < \phi < 1/2:
    I = (4 \pi t / N) \sin(2 \pi (1/2 - \phi) / N) \sin(2 \pi m / N) /
         \sin(\pi / N?)
 - Two-terminal transmission: semi-infinite tight-binding leads (hopping t)
   attached at sites a and b
   With \Sigma(E) = (E - i \sqrt(4 t^2 - E^2)) / 2 and
        \Gamma = \sqrt(4 t^2 - E^2)
   Green's function G = (E - H - \Sigma_a - \Sigma_b)^-1
   gives T = \Gamma^2 |G_ab|^2 and R = |1 - i \Gamma G_aa|^2 and T + R = 1
   (unitarity)
   2 arms interfere with a phase difference (2 \pi \phi) (+ dynamical phase
   of arm-length difference), so T oscillates with period one flux quantum
```

## Types

### `ar_scatter_t`

```c
typedef struct {
  double transmission; /* |S_ba|^2 */
  double reflection;   /* |S_aa|^2 */
} ar_scatter_t;
```

## Functions (5)

[`ar_level`](#ar_level), [`ar_levels_sorted`](#ar_levels_sorted), [`ar_energy`](#ar_energy), [`ar_current`](#ar_current), [`ar_scatter`](#ar_scatter)

### `ar_level`

```c
double ar_level(int n_sites, int n, double flux, double t);
```

```text
E_n of level n (any integer n; periodic)
```

### `ar_levels_sorted`

```c
int ar_levels_sorted(int n_sites, double flux, double t, double *e);
```

```text
All N levels in ascending order (stable on ties). e must hold n_sites
doubles; returns 0 or -1 on bad input
```

### `ar_energy`

```c
double ar_energy(int n_sites, int ne, double flux, double t);
```

```text
Total energy of ne lowest levels (0 <= ne <= N); NaN on bad input
```

### `ar_current`

```c
double ar_current(int n_sites, int ne, double flux, double t);
```

```text
Persistent current -dE/dphi of ne lowest levels (analytic derivative); NaN on
bad input
```

### `ar_scatter`

```c
int ar_scatter(int n_sites, int a, int b, double flux, double energy, double t,
               ar_scatter_t *out);
```

```text
Transmission and reflection for leads at sites a != b at energy |E| < 2t; -1
on bad input. 4 <= N <= 64
```
