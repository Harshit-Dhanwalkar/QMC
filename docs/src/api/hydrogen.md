# `hydrogen`

> Generated from `physics/hydrogen.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Hydrogen atom: radial Schrodinger solver, energy levels and radial wavefunctions.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/hydrogen.h>
```

**Start with:** [`hydrogen_radial_solve`](#hydrogen_radial_solve), [`hydrogen_energy_level`](#hydrogen_energy_level), [`hydrogen_radial_wavefunction`](#hydrogen_radial_wavefunction)

- Guide: [physics/hydrogen.md](../physics/hydrogen.md)
- Example: `examples/eg_06_hydrogen.c`
- Test (every call with expected results): `tests/test_hydrogen.c`
- Runtime help: `qmc_help(stdout, "hydrogen");` or `qmc-help hydrogen`

## Functions (3)

[`hydrogen_radial_solve`](#hydrogen_radial_solve), [`hydrogen_energy_level`](#hydrogen_energy_level), [`hydrogen_radial_wavefunction`](#hydrogen_radial_wavefunction)

### `hydrogen_radial_solve`

```c
eigen_t *hydrogen_radial_solve(const double *r, int N, int l, double hbar,
                               double mass, double e_charge, double eps0);
```

```text
Radial Schrodinger equation for hydrogen:
[ -\hbar^2 / (2m) * (d^2 / dr^2 - l(l+1) / r^2) - e^2 / (4 * \pi * \epsilon_0
r) ] R(r) = E R(r)

Returns eigen_t with eigenvalues (energies) and eigenvectors (radial
functions).
NOTE: The radial grid r[0..N-1] must be given (typically logarithmic). l is
the angular momentum quantum number.
```

### `hydrogen_energy_level`

```c
double hydrogen_energy_level(int n);
```

```text
Analytic hydrogen energy levels (SI units)
```

### `hydrogen_radial_wavefunction`

```c
cvector_t *hydrogen_radial_wavefunction(const double *r, int N, int n, int l);
```

```text
Hydrogen radial wavefunction (analytic, unnormalized) for given n, l.

Returns new cvector_t* with values on the grid r.
Where
- r: grid
- N: size
- n: principal quantum number
- l: angular momentum
NOTE: The function uses associated Laguerre polynomials.
```
