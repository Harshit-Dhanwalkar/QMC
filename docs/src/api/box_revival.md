# `box_revival`

> Generated from `physics/box_revival.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Particle in a box: Gaussian packet evolution, wavefunction revivals and density maps.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/box_revival.h>
```

**Start with:** [`box_gaussian_coefficients`](#box_gaussian_coefficients), [`box_wavefunction`](#box_wavefunction), [`box_density_map`](#box_density_map), [`box_overlap`](#box_overlap), [`box_mean_energy`](#box_mean_energy)

- Example: `examples/eg_89_box_revival.c`
- Test (every call with expected results): `tests/test_box_revival.c`
- Runtime help: `qmc_help(stdout, "box_revival");` or `qmc-help box_revival`

## Overview

```text
Quantum revivals and carpets of a wave packet in an infinite square well

Units: well [0, 1], \hbar = 1, mass \pi^2 / 2 so that E_n = n^2 and
eigenfunctions are \phi_n(x) = \sqrt{2} \sin(n \pi x)

A state \sum_n c_n \phi_n evolves as c_n \exp(-i n^2 t), so it is periodic
with revival time T = 2 \pi: at t = T packet is exactly back, at t = T/2 it
is minus mirror image about x = 1/2, \psi(x, T/2) = -\psi_0(1 - x) and at t =
T/4 it is superposition ((1 - i) / 2) \psi_0(x) - ((1 + i) / 2) \psi_0(1 - x)
of packet and its mirror image (a cat), and so on for every rational fraction
of T (Talbot carpet)
```

## Macros and constants

```c
#define BOX_MAX_MODES 512
```

```c
#define BOX_REVIVAL_TIME 6.283185307179586
```

## Functions (6)

[`box_gaussian_coefficients`](#box_gaussian_coefficients), [`box_wavefunction`](#box_wavefunction), [`box_density_map`](#box_density_map), [`box_overlap`](#box_overlap), [`box_mean_energy`](#box_mean_energy), [`box_mean_position`](#box_mean_position)

### `box_gaussian_coefficients`

```c
double box_gaussian_coefficients(complex_t *c, int n_modes, double x0,
                                 double sigma, double n0);
```

```text
Project a Gaussian packet onto first n_modes eigenstates:
 \psi(x) ~ \exp(-(x - x0)^2 / (4 \sigma^2) + i k0 x)
 \sigma = position spread
 k0 = \pi * n0
so n0 is central mode number. Coefficients are then renormalised to unit
norm. Returns squared norm that fell inside first n_modes before
renormalisation (1 means no truncation error), or a negative number on bad
input (0 < x0 < 1, 0.005 <= sigma <= 0.25, 1 <= n_modes <= BOX_MAX_MODES)
```

### `box_wavefunction`

```c
int box_wavefunction(const complex_t *c, int n_modes, double t, const double *x,
                     int nx, complex_t *psi);
```

```text
Wavefunction \psi(x, t) = \sum c_n \exp(-i n^2 t) \phi_n(x) at nx points. 0
on success, -1 on bad input
```

### `box_density_map`

```c
int box_density_map(const complex_t *c, int n_modes, const double *x, int nx,
                    double t_max, int nt, double *rho);
```

```text
Space-time density |\psi(x_j, t_i)|^2 for t_i = t_max i / (nt - 1), written
to \rho[i * nx + j]. Same physics as box_wavefunction but with sine table
computed once and negligible modes (|c_n|^2 < 1e-16) skipped. 0 on success,
-1 on bad input or allocation failure
```

### `box_overlap`

```c
complex_t box_overlap(const complex_t *c, int n_modes, double t);
```

```text
Return amplitude <\psi(0)|\psi(t)> = \sum |c_n|^2 \exp(-i n^2 t)
```

### `box_mean_energy`

```c
double box_mean_energy(const complex_t *c, int n_modes);
```

```text
Mean energy \sum n^2 |c_n|^2
```

### `box_mean_position`

```c
double box_mean_position(const double *x, const complex_t *psi, int nx);
```

```text
Mean position <x> from wavefunction on grid (trapezoid rule)
```
