# `numerov`

> Generated from `core/ode/numerov.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Numerov integrator for the 1D time-independent Schrodinger equation.

**Category:** Core: ODE & time evolution

```c
#include <qmc/core/ode/numerov.h>
```

**Start with:** [`numerov_shoot`](#numerov_shoot), [`numerov_integrate`](#numerov_integrate), [`numerov_shoot_matching`](#numerov_shoot_matching)

- Guide: [internals/numerov.md](../internals/numerov.md)
- Test (every call with expected results): `tests/test_numerov.c`
- Runtime help: `qmc_help(stdout, "numerov");` or `qmc-help numerov`

## Overview

```text
Numerov algorithm for 1D time-independent Schr$\"{o}$dinger equation:
  -\hbar^2/(2m) d^2\phi/dx^2 + V(x)\phi = E \cdot \phi
```

## Types

### `numerov_params_t`

```c
typedef struct {
  double *x;         // Position grid
  double *V;         // Potential array V(x)
  int n;             // Grid points
  double dx;         // Grid spacing
  double hbar_sq_2m; // \hbar^2/(2m) in problem units
} numerov_params_t;
```

### `numerov_solution_t`

```c
typedef struct {
  double energy;
  cvector_t *psi;
} numerov_solution_t;
```

## Functions (4)

[`numerov_shoot`](#numerov_shoot), [`numerov_integrate`](#numerov_integrate), [`numerov_solution_free`](#numerov_solution_free), [`numerov_shoot_matching`](#numerov_shoot_matching)

### `numerov_shoot`

```c
numerov_solution_t *numerov_shoot(numerov_params_t *params, double E_guess,
                                  double E_tol);
```

```text
numerov_shoot: find eigenstate for given level via matrix
diagonalization (tridiag_eigh) of discretized Hamiltonian.

E_guess is used only to determine target level index:
  level = round(E_guess - V_min - 0.5), clamped to >= 0.
```

### `numerov_integrate`

```c
void numerov_integrate(const numerov_params_t *params, double E,
                       cvector_t *psi);
```

```text
numerov_integrate: Numerov integration at given energy.

Fills \psi with wavefunction by forward Numerov integration.
Seed: \psi[0]=0 (Dirichlet), psi[1]=1e-8.
```

### `numerov_solution_free`

```c
void numerov_solution_free(numerov_solution_t *sol);
```

_No description in the header yet._

### `numerov_shoot_matching`

```c
numerov_solution_t *numerov_shoot_matching(numerov_params_t *params,
                                           double E_min, double E_max,
                                           int n_scan, double tol);
```

```text
numerov_shoot_matching: find a bound-state eigenvalue via Numerov shooting
with log-derivative matching at outer classical turning point,

params:       same grid/potential convention as numerov_shoot()
E_min, E_max: energy bracket containing exactly one eigenvalue
n_scan:       number of points used to scan [E_min,E_max] for a smooth
              sign change before bisecting (200 is a reasonable default)
tol:          energy convergence tolerance for bisection
```
