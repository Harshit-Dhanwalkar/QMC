# `bloch`

> Generated from `physics/bloch.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Single-qubit Bloch-vector picture on top of the Lindblad solver.

**Category:** Two-level, open systems & quantum optics

```c
#include <qmc/physics/bloch.h>
```

**Start with:** [`bloch_vector_from_density`](#bloch_vector_from_density), [`bloch_density_from_vector`](#bloch_density_from_vector), [`bloch_evolve`](#bloch_evolve)

- Example: `examples/eg_78_bloch.c`
- Test (every call with expected results): `tests/test_bloch.c`
- Runtime help: `qmc_help(stdout, "bloch");` or `qmc-help bloch`

## Overview

```text
Single-qubit Bloch-vector picture built on Lindblad solver (lindblad.h)

Basis: index 0 = ground state |g> (north pole, z = +1), index 1 = excited
state |e> (south pole, z = -1)
Density matrix is \rho = (1 + v . \sigma) / 2 with
   x = 2 Re \rho_01,   y = -2 Im \rho_01,   z = \rho_00 - \rho_11,
|v| = 1 for pure states and |v| < 1 for mixed ones

Driven, damped qubit (\hbar = 1), RWA Hamiltonian as in rabi.h
    H = (1/2) [[ \Delta, \Omega ], [ \Omega, -\Delta ]] = (1/2) b . \sigma,
    b = (\Omega, 0, \Delta),
so without damping v precesses about b at generalised Rabi frequency
|b| = \sqrt(\Omega^2 + \Delta^2): dv/dt = b x v

Decoherence enters through jump operators from lindblad.h:
  gamma1    amplitude damping L = sqrt(gamma1) |g><e|     (T1 = 1/gamma1)
  gamma_phi pure dephasing    L = sqrt(gamma_phi/2) sigma_z
giving, with no drive, z -> 1 - (1 - z0) e^{-gamma1 t} and transverse
components decaying at rate gamma1/2 + gamma_phi (1/T2).
```

## Functions (3)

[`bloch_vector_from_density`](#bloch_vector_from_density), [`bloch_density_from_vector`](#bloch_density_from_vector), [`bloch_evolve`](#bloch_evolve)

### `bloch_vector_from_density`

```c
int bloch_vector_from_density(const cmatrix_t *rho, double v[3]);
```

```text
Bloch vector of a 2x2 density matrix. Returns 0, or -1 if rho is NULL or
not 2x2.
```

### `bloch_density_from_vector`

```c
cmatrix_t *bloch_density_from_vector(const double v[3]);
```

```text
New 2x2 density matrix (1 + v . sigma)/2 for |v| <= 1 (small round-off
tolerance); NULL for non-finite v, |v| > 1 + 1e-9, or allocation failure.
Caller frees with cmatrix_free.
```

### `bloch_evolve`

```c
int bloch_evolve(double v[3], double omega, double delta, double gamma1,
                 double gamma_phi, double dt, int steps);
```

```text
Advance Bloch vector v (in place) by `steps` RK4 steps of size dt under drive
\omega, detuning \delta and damping rates above

Returns 0 on success, -1 on invalid input (NULL v, non-finite entries,
|v| > 1 + 1e-9, non-finite \omega/\delta, negative or non-finite rate,
dt <= 0, steps < 1), -2 on allocation failure. On failure v is untouched
```
