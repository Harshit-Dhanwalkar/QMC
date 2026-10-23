# `geometry_optimization`

> Generated from `physics/geometry_optimization.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Generic Cartesian geometry optimiser driven by an energy-and-gradient callback.

**Category:** Quantum chemistry (molecules)

```c
#include <qmc/physics/geometry_optimization.h>
```

**Start with:** [`optimize_geometry_steepest_descent`](#optimize_geometry_steepest_descent)

- Example: `examples/eg_49_geometry_optimization.c`
- Test (every call with expected results): `tests/test_geometry_optimization.c`
- Runtime help: `qmc_help(stdout, "geometry_optimization");` or `qmc-help geometry_optimization`

## Overview

```text
Generic Cartesian geometry optimizer, driven by an energy+analytic-gradient
callback (Not tied to RHF/UHF or to any one molecule's basis-set
construction)
NOTE: eg_49_geometry_optimization.c and eg_66_uhf_geometry_optimization.c
each hand-roll a scalar (bond-length-only) steepest-descent loop that
rebuilds whole basis/molecule per force evaluation; this generalizes that
idea to a full n_atoms*3 Cartesian coordinate vector so it also covers
non-collinear (bent, polyatomic) geometries, while keeping same "rebuild
from these coordinates and hand back energy+gradient" shape those examples
already use
```

## Types

### `geom_energy_grad_fn`

```text
Evaluate energy and its Cartesian gradient at a given geometry
 coords    : input, length n_atoms*3 (coords[3 * A + d], d=0,1,2 for x,y,z)
 n_atoms   : number of atoms (fixed for lifetime of optimization)
 user_data : opaque pointer forwarded from optimizer call, typically
             holding basis-set kind, charges, electron count, and any SCF
             settings needed to rebuild basis at new geometry
energy_out : callee writes total (electronic + nuclear repulsion) energy
             at coords
grad_out   : preallocated length n_atoms*3 by caller; callee writes
             dE/dR_A[d] into grad_out[3*A+d]

Returns 0 on success (e.g. SCF converged), nonzero on failure (e.g. SCF did
not converge at this geometry) - optimizer aborts on a nonzero return and
reports it as not converged
```

```c
typedef int (*geom_energy_grad_fn)(const double *coords, int n_atoms,
                                   void *user_data, double *energy_out,
                                   double *grad_out);
```

### `geometry_optimization_result_t`

```c
typedef struct {
  double *coords;   /* final geometry, length n_atoms*3, caller frees */
  double energy;    /* energy at final geometry */
  double grad_norm; /* max |component| of gradient at final geometry */
  int iterations;
  int converged;
} geometry_optimization_result_t;
```

## Functions (2)

[`optimize_geometry_steepest_descent`](#optimize_geometry_steepest_descent), [`geometry_optimization_result_free`](#geometry_optimization_result_free)

### `optimize_geometry_steepest_descent`

```c
geometry_optimization_result_t *
optimize_geometry_steepest_descent(geom_energy_grad_fn func, void *user_data,
                                   const double *coords0, int n_atoms,
                                   double step0, double grad_tol, int max_iter);
```

```text
Steepest descent on analytic gradient, with simple backtracking (halving
step) whenever a step would raise energy - same "move downhill
along force", made robust to a poorly-scaled initial step and generalized
to n_atoms*3 coordinates instead of one scalar bond length

func/user_data: see geom_energy_grad_fn above
coords0       : starting geometry, length n_atoms*3 (not modified)
n_atoms       : number of atoms
step0         : initial step length along -gradient
grad_tol      : converged when every |grad component| < grad_tol
max_iter      : outer iteration cap (each outer iteration may backtrack
                step internally without consuming an outer iteration)

Returns NULL on invalid input (n_atoms <= 0, step0 <= 0, grad_tol <= 0,
max_iter <= 0) or allocation failure. Otherwise always returns a result
(check -> converged)
```

### `geometry_optimization_result_free`

```c
void geometry_optimization_result_free(geometry_optimization_result_t *res);
```

_No description in the header yet._
