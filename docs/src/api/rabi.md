# `rabi`

> Generated from `physics/rabi.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Rabi oscillations of a two-level system in the rotating-wave approximation, exact and numerical.

**Category:** Two-level, open systems & quantum optics

```c
#include <qmc/physics/rabi.h>
```

**Start with:** [`rabi_excited_probability`](#rabi_excited_probability), [`rabi_evolve_exact`](#rabi_evolve_exact), [`rabi_evolve_numerical`](#rabi_evolve_numerical)

- Guide: [physics/rabi.md](../physics/rabi.md)
- Example: `examples/eg_14_rabi.c`
- Test (every call with expected results): `tests/test_rabi.c`
- Runtime help: `qmc_help(stdout, "rabi");` or `qmc-help rabi`

## Overview

```text
Two-level system in the rotating-wave approximation (RWA).

H = (\hbar/2) * [[ \Delta, \Omega ], [ \Omega, -\Delta ]]
in the {|1> (ground), |2> (excited)} basis. \Delta = \omega_drive -
\omega_0 is the detuning (rad/s) of drive from bare transition
frequency; \Omega is the real, angular-frequency-units Rabi coupling
strength.

H is time-independent in rotating frame
```

## Functions (3)

[`rabi_excited_probability`](#rabi_excited_probability), [`rabi_evolve_exact`](#rabi_evolve_exact), [`rabi_evolve_numerical`](#rabi_evolve_numerical)

### `rabi_excited_probability`

```c
double rabi_excited_probability(double t, double Omega, double Delta);
```

```text
Exact excited-state probability starting from ground state at
t=0: P_e(t) = (\Omega^2 / \Omega_R^2) * \sin^2(\Omega_R*t/2), where
\Omega_R = \sqrt(\Omega^2 + \Delta^2) is generalized Rabi frequency.
```

### `rabi_evolve_exact`

```c
int rabi_evolve_exact(cvector_t *psi, double t, double Omega, double Delta);
```

```text
Exact unitary propagation of two-level spinor \psi (length 2, \psi[0]=ground
amplitude, \psi[1]=excited amplitude) from t=0 to time t, using closed-form
evolution operator

U(t) = \cos(\Omega_R t/2) I - i * \sin(\Omega_R * t/2)/\Omega_R *
[[\Delta,\Omega],[\Omega ,-\Delta]]

Overwrites psi in place.

Returns 0 on success, -1 on invalid input.
```

### `rabi_evolve_numerical`

```c
int rabi_evolve_numerical(cvector_t *psi, double hbar, double Omega,
                          double Delta, double dt, int steps);
```

```text
Treating H as a genuine 2x2 real symmetric (trivially tridiagonal) energy
Hamiltonian: diag = [\hbar * \Delta/2, -\hbar * \Delta/2], offdiag = [\hbar *
\Omega/2]
```
