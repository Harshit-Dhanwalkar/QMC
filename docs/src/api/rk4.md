# `rk4`

> Generated from `core/ode/rk4.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Fourth-order Runge-Kutta integrator for complex ODE systems.

**Category:** Core: ODE & time evolution

```c
#include <qmc/core/ode/rk4.h>
```

**Start with:** [`rk4_step`](#rk4_step)

- Test (every call with expected results): `tests/test_rk4.c`
- Runtime help: `qmc_help(stdout, "rk4");` or `qmc-help rk4`

## Overview

```text
Runge-Kutta 4th order integrator for complex ODEs.
System: d\phi/dt = f(\phi, t).
```

## Types

### `ode_func_t`

```c
typedef void (*ode_func_t)(double t, const cvector_t *y, cvector_t *dydt,
                           void *params);
```

## Functions (1)

[`rk4_step`](#rk4_step)

### `rk4_step`

```c
int rk4_step(double t, double dt, cvector_t *y, ode_func_t f, void *params);
```

```text
Perform one RK4 step from time t to t+dt.
y: input/output vector (complex)
params: additional parameters for f.
Returns 0 on success.
```
