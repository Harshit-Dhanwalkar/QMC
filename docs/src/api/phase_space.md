# `phase_space`

> Generated from `physics/phase_space.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Phase-space picture of a bosonic mode: Fock states, Kerr nonlinearity and the Wigner function.

**Category:** Two-level, open systems & quantum optics

```c
#include <qmc/physics/phase_space.h>
```

**Start with:** [`ps_fock_state`](#ps_fock_state), [`ps_coherent_state`](#ps_coherent_state), [`ps_cat_state`](#ps_cat_state), [`ps_squeezed_vacuum`](#ps_squeezed_vacuum), [`ps_kerr_evolve`](#ps_kerr_evolve)

- Example: `examples/eg_86_phase_space.c`
- Test (every call with expected results): `tests/test_phase_space.c`
- Runtime help: `qmc_help(stdout, "phase_space");` or `qmc-help phase_space`

## Overview

```text
Phase-space picture of a single bosonic mode: Fock-space states, Kerr
nonlinearity and Wigner function

A pure state is stored as amplitudes c[0 .. n_max] in Fock basis,
|\psi> = \sum_n c_n |n>, truncated at n_max. Units: \hbar = m = \omega = 1,
so x = (a + a^\dagger) / \sqrt(2), p = (a - a^\dagger) / (i \sqrt(2)) and
coherent state |\alpha> sits at (x, p) = \sqrt(2) (Re \alpha, Im \alpha) with
a Gaussian Wigner function exp(-(x - x0)^2 - (p - p0)^2) / \pi

Kerr Hamiltonian H = \chi n (n - 1) is diagonal in Fock basis, so its
evolution is exact: c_n -> c_n \exp(-i \chi t n (n - 1)). It twists phase
space by an amount that grows with n, so a coherent state first smears into a
banana, at \chi t = \pi / 2 splits into a Schroedinger-cat superposition of
|i \alpha> and |-i \alpha> (exactly, with relative phase \exp(-i \pi / 2) and
at \chi t = \pi returns to itself (full revival)

Wigner function is evaluated from position wave function
  W(x, p) = (1 / \pi) \int \psi*(x + y) \psi(x - y) \exp(2 i p y) dy
on a uniform x grid, with integral done as a sum over grid offsets
It is negative wherever state has non-classical interference:
Fock states n >= 1, cat states, squeezed states (no negativity but sub-vacuum
variance)
```

## Functions (10)

[`ps_fock_state`](#ps_fock_state), [`ps_coherent_state`](#ps_coherent_state), [`ps_cat_state`](#ps_cat_state), [`ps_squeezed_vacuum`](#ps_squeezed_vacuum), [`ps_kerr_evolve`](#ps_kerr_evolve), [`ps_mean_photon_number`](#ps_mean_photon_number), [`ps_norm2`](#ps_norm2), [`ps_wavefunction`](#ps_wavefunction), [`ps_wigner`](#ps_wigner), [`ps_wigner_negativity`](#ps_wigner_negativity)

### `ps_fock_state`

```c
int ps_fock_state(complex_t *c, int n_max, int n);
```

```text
Fock state |n>

Returns 0, or -1 for NULL c, n_max < 0 or n outside [0, n_max]
```

### `ps_coherent_state`

```c
int ps_coherent_state(complex_t *c, int n_max, double alpha_re,
                      double alpha_im);
```

```text
Coherent state |\alpha>, normalised on truncated space

Returns 0, or -1 for NULL c, n_max < 0, non-finite alpha or |\alpha|^2 above
n_max / 2 (truncation would no longer hold state)
```

### `ps_cat_state`

```c
int ps_cat_state(complex_t *c, int n_max, double alpha_re, double alpha_im,
                 int parity);
```

```text
NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
```

### `ps_squeezed_vacuum`

```c
int ps_squeezed_vacuum(complex_t *c, int n_max, double r, double theta);
```

```text
Squeezed vacuum S(r \exp(i \theta)) |0>: <x^2> = \exp(-2r) / 2 for \theta = 0

Returns 0, or -1 for NULL c, n_max < 0, non-finite r or \theta, or a state
whose mean photon number \sinh^2(r) exceeds n_max / 3
```

### `ps_kerr_evolve`

```c
int ps_kerr_evolve(complex_t *c, int n_max, double chi_t);
```

```text
Evolve under H = \chi n (n - 1) for \chi t = \chi_t

Returns 0, or -1 for NULL c, n_max < 0 or non-finite \chi_t
```

### `ps_mean_photon_number`

```c
double ps_mean_photon_number(const complex_t *c, int n_max);
```

```text
Mean photon number and squared norm. NaN for NULL c or n_max < 0
```

### `ps_norm2`

```c
double ps_norm2(const complex_t *c, int n_max);
```

_No description in the header yet._

### `ps_wavefunction`

```c
int ps_wavefunction(const complex_t *c, int n_max, const double *x, int nx,
                    complex_t *psi);
```

```text
Position wave function \psi(x_j) = \sum_n c_n \phi_n(x_j) at nx points, with
Hermite functions phi_n from stable three-term recurrence

Returns 0, or -1 for NULL arguments, n_max < 0 or nx < 1
```

### `ps_wigner`

```c
int ps_wigner(const complex_t *psi, int nx, double dx, const double *p, int np,
              double *W);
```

```text
Wigner function on grid x_j = x0 + j dx (j = 0 .. nx - 1) and nominal momenta
p[0 .. np - 1]: W[l * nx + j] = W(x_j, p_l). \psi holds wave function on x
grid. p grid should stay within \pi / (2 dx) to avoid aliasing

Returns 0, or -1 for NULL arguments, nx < 1, np < 1 or dx <= 0
```

### `ps_wigner_negativity`

```c
double ps_wigner_negativity(const double *W, int nx, int np, double dx,
                            double dp);
```

```text
Negativity (integral |W| - 1) / 2 of a Wigner grid with spacings dx, dp:
total volume of negative part. NaN for NULL W or non-positive sizes
```
