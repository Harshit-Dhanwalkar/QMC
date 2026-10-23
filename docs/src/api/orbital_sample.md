# `orbital_sample`

> Generated from `physics/orbital_sample.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Monte Carlo sampling of hydrogen orbitals psi_nlm, for visualisation.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/orbital_sample.h>
```

**Start with:** [`hydrogen_orbital_sample`](#hydrogen_orbital_sample), [`hydrogen_orbital_radial_au`](#hydrogen_orbital_radial_au), [`hydrogen_orbital_angular_density`](#hydrogen_orbital_angular_density)

- Example: `examples/eg_77_orbital_sample.c`
- Test (every call with expected results): `tests/test_orbital_sample.c`
- Runtime help: `qmc_help(stdout, "orbital_sample");` or `qmc-help orbital_sample`

## Overview

```text
Monte Carlo sampling of hydrogen orbitals \psi_{nlm} = R_nl(r) Y_lm(\theta,
\phi)

Positions are drawn from |\psi|^2 exactly (up to table resolution), so a
scatter plot of points is orbital: density of dots = probability density.
All lengths here are in Bohr radii a0 and R_nl is in a0^(-3/2)

Radial part comes from library's analytic hydrogen_radial_wavefunction()
(hydrogen.h); radial distribution r^2 R^2 is tabulated and sampled by
inverse CDF. Angular part uses library's spherical harmonics and rejection
sampling in (\cos(\theta, \phi)) (uniform in \cos(\theta) carries
\sin(\theta) volume element)

real_form = 0: complex orbitals Y_l^m (density independent of \phi)
real_form = 1: real combinations:
    m = 0: Y_l0,  m > 0: \sqrt(2) Re Y_l^m (\cos(m \phi)),
    m < 0: \sqrt(2) Im Y_l^|m| (\sin |m| \phi)
```

## Types

### `?`

```c
enum { ORBITAL_MAX_N = 8 };
```

## Functions (3)

[`hydrogen_orbital_sample`](#hydrogen_orbital_sample), [`hydrogen_orbital_radial_au`](#hydrogen_orbital_radial_au), [`hydrogen_orbital_angular_density`](#hydrogen_orbital_angular_density)

### `hydrogen_orbital_sample`

```c
int hydrogen_orbital_sample(int n, int l, int m, int real_form, int count,
                            unsigned long long seed, double *out);
```

```text
Draw `count` points from |\psi_{nlm}|^2. `out` must hold 4*count doubles;
i-th point is written as (x, y, z, phase) at out[4i .. 4i+3]:
  x, y, z : position in a0
  phase   : complex form: arg \psi in (-\pi, \pi]; real form: 0 where \psi >=
            0, \pi where \psi < 0 (the sign lobes of a real orbital)

The sequence is deterministic for a given `seed` (own xorshift generator)

Returns 0 on success, -1 on invalid input (n outside 1..ORBITAL_MAX_N, l
outside 0..n-1, |m| > l, count < 1, out NULL), -2 on allocation failure
```

### `hydrogen_orbital_radial_au`

```c
double hydrogen_orbital_radial_au(int n, int l, double r);
```

```text
Radial function R_nl(r) in a0^(-3/2) at r (in a0), normalised so that \int
R^2 r^2 dr = 1. NaN for invalid n, l or r < 0
```

### `hydrogen_orbital_angular_density`

```c
double hydrogen_orbital_angular_density(int l, int m, int real_form,
                                        double theta, double phi);
```

```text
Angular probability density |Y(\theta, \phi)|^2 per steradian of complex
(real_form = 0) or real (real_form = 1) harmonic. NaN for invalid l, m
```
