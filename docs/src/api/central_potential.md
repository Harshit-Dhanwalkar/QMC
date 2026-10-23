# `central_potential`

> Generated from `physics/central_potential.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

General 3D central-potential radial solver.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/central_potential.h>
```

**Start with:** [`central_potential_radial_solve`](#central_potential_radial_solve)

- Guide: [physics/central_potential.md](../physics/central_potential.md)
- Example: `examples/eg_07_central_potential.c`
- Test (every call with expected results): `tests/test_central_potential.c`
- Runtime help: `qmc_help(stdout, "central_potential");` or `qmc-help central_potential`

## Functions (1)

[`central_potential_radial_solve`](#central_potential_radial_solve)

### `central_potential_radial_solve`

```c
eigen_t *central_potential_radial_solve(const double *r, int N, int l,
                                        double hbar, double mass,
                                        potential_fn V, void *params);
```

```text
General 3D central-potential radial solver.

Solves the radial Schrodinger equation for u(r) = r * R(r):
  -\hbar^2/(2m) u''(r) + [ V(r) + \hbar^2/(2m) * l(l+1)/r^2 ] u(r) = E u(r)

on a uniform grid r[0..N-1], via finite-difference tridiagonal
discretization + tridiag_eigh (implicit QL).
```
