# `double_well`

> Generated from `physics/double_well.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Symmetric double-well tunnelling: parity eigenstates, exponential splitting versus barrier, coherent oscillation between wells.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/double_well.h>
```

**Start with:** [`dw_potential`](#dw_potential), [`dw_eigen`](#dw_eigen), [`dw_project`](#dw_project), [`dw_evolve`](#dw_evolve)

- Example: `examples/eg_96_double_well.c`
- Test (every call with expected results): `tests/test_double_well.c`
- Runtime help: `qmc_help(stdout, "double_well");` or `qmc-help double_well`

## Overview

```text
Tunnelling in a symmetric double well (ammonia-inversion physics)

 H = -(1/2) d^2/dx^2 + V(x)
with \hbar = m = 1 on a
 Dirichlet grid x_i = -L + (i + 1) dx
 dx = 2L / (n + 1)
 V(x) = B (1 - x^2 / a^2)^2

two wells at x = -+a separated by a barrier of height B

Lowest eigenstates come in parity doublets (even, odd) whose splitting
\Delta = E1 - E0 shrinks exponentially with barrier
(WKB: \Delta ~ (\omega / \pi) \exp(-\int \sqrt(2 (V - E)) dx across barrier,
\omega = \sqrt(8 B / a^2) small-oscillation frequency of each well)
A particle prepared in one well, (\psi_0 - \psi_1) / \sqrt(2) (with sign
convention below), tunnels to other and back with P_left(t) = \cos^2(\Delta t
/ 2): it oscillates with period 2 \pi / \Delta

Solver is a tridiagonal Sturm bisection with parity-projected inverse
iteration, so near-degenerate doublets come out cleanly; potential must be
even
```

## Types

### `?`

```c
enum { DW_MAX_GRID = 4096, DW_MAX_STATES = 12 };
```

## Functions (7)

[`dw_potential`](#dw_potential), [`dw_spacing`](#dw_spacing), [`dw_eigen`](#dw_eigen), [`dw_project`](#dw_project), [`dw_evolve`](#dw_evolve), [`dw_prob_left`](#dw_prob_left), [`dw_mean_x`](#dw_mean_x)

### `dw_potential`

```c
int dw_potential(double barrier, double a, double half_width, int n, double *v);
```

```text
Fills v with V(x_i) = B (1 - x^2 / a^2)^2 on grid of n points over
[-half_width, half_width]; 0, or -1 on bad input
```

### `dw_spacing`

```c
double dw_spacing(int n, double half_width);
```

```text
Grid spacing for n points on [-L, L].
```

### `dw_eigen`

```c
int dw_eigen(const double *v, int n, double dx, int n_states, double *e,
             double *psi);
```

```text
Lowest n_states eigenvalues e[k] (ascending) and normalised eigenvectors
\psi[k * n + i] (\sum \psi^2 dx = 1) of an even potential v. Even states are
positive at their central-most maximum; odd states are positive on right
side. 0 on success, -1 on bad input (n <= DW_MAX_GRID,
n_states <= DW_MAX_STATES, v even to 1e-9 relative)
```

### `dw_project`

```c
int dw_project(const double *psi, int n, double dx, int n_states,
               double half_width, double center, double sigma, double *c,
               double *captured);
```

```text
Overlaps c[k] = <\psi_k | g> of a normalised Gaussian g(x) ~ \exp(-(x -
center)^2 / (4 \sigma^2)) with eigenstates; *captured receives sum c^2 (1
when states span packet). 0 or -1
```

### `dw_evolve`

```c
int dw_evolve(const double *psi, const double *e, const double *c, int n,
              int n_states, double t, double *re, double *im);
```

```text
\psi(t) = \sum_k c_k \exp(-i e_k t) \psi_k as (re, im) arrays of length n. 0
or -1
```

### `dw_prob_left`

```c
double dw_prob_left(const double *re, const double *im, int n, double dx);
```

```text
Probability to left of barrier (x < 0) of a state (re, im), normalised by
total probability
```

### `dw_mean_x`

```c
double dw_mean_x(const double *re, const double *im, int n, double dx,
                 double half_width);
```

```text
Mean position <x> for a state (re, im), normalised
```
