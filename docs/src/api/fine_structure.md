# `fine_structure`

> Generated from `physics/fine_structure.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Hydrogen fine structure: spin-orbit coupling and the closed-form total fine-structure formula.

**Category:** Atomic structure & perturbation theory

```c
#include <qmc/physics/fine_structure.h>
```

**Start with:** [`spin_orbit_ls_expect`](#spin_orbit_ls_expect), [`hydrogen_expect_inv_r3`](#hydrogen_expect_inv_r3), [`hydrogen_spin_orbit_energy`](#hydrogen_spin_orbit_energy), [`hydrogen_fine_structure_shift`](#hydrogen_fine_structure_shift), [`spin_orbit_ls_expect_from_coupling`](#spin_orbit_ls_expect_from_coupling)

- Guide: [physics/fine_structure.md](../physics/fine_structure.md)
- Example: `examples/eg_16_finestructure.c`
- Test (every call with expected results): `tests/test_fine_structure.c`
- Runtime help: `qmc_help(stdout, "fine_structure");` or `qmc-help fine_structure`

## Overview

```text
Hydrogen fine structure: spin-orbit coupling closed-form total fine-structure
formula for cross-checking.

Conventions: l is the orbital quantum number (integer, l >= 0). Spin is fixed
at s=1/2 (electron).
```

## Functions (5)

[`spin_orbit_ls_expect`](#spin_orbit_ls_expect), [`hydrogen_expect_inv_r3`](#hydrogen_expect_inv_r3), [`hydrogen_spin_orbit_energy`](#hydrogen_spin_orbit_energy), [`hydrogen_fine_structure_shift`](#hydrogen_fine_structure_shift), [`spin_orbit_ls_expect_from_coupling`](#spin_orbit_ls_expect_from_coupling)

### `spin_orbit_ls_expect`

```c
double spin_orbit_ls_expect(int l, int j_2);
```

```text
<L.S>/\hbar^2 = (1/2) * (j(j+1) - l(l+1) - s(s+1)), s=1/2 fixed. Closed form
```

### `hydrogen_expect_inv_r3`

```c
double hydrogen_expect_inv_r3(int n, int l, double hbar, double mass,
                              double e_charge, double eps0);
```

```text
<1/r^3>_{nl} for hydrogen = 1/(a0^3 n^3 l(l+1/2)(l+1)).
Returns 0 for l=0
```

### `hydrogen_spin_orbit_energy`

```c
double hydrogen_spin_orbit_energy(int n, int l, int j_2, double hbar,
                                  double mass, double e_charge, double eps0,
                                  double c);
```

```text
Perturbative spin-orbit energy shift for hydrogen:
  dE_SO = \exp^2/(8 * \pi * \eps0 m^2 * c^2) * <1/r^3>_{nl} * \hbar^2 *
          <L.S>/\hbar^2

Returns 0 for l=0.

NOTE: spin-orbit contribution only. Total physical fine-structure shift also
includes relativistic kinetic-energy correction and (for l=0) Darwin term
```

### `hydrogen_fine_structure_shift`

```c
double hydrogen_fine_structure_shift(int n, int j_2, double hbar, double mass,
                                     double e_charge, double eps0, double c);
```

```text
Known closed-form total fine-structure shift:
  dE_fs = E_n * (\alpha^2 / n^2) * (n / (j+1/2) - 3/4)
Where E_n = hydrogen_energy_level(n)
```

### `spin_orbit_ls_expect_from_coupling`

```c
double spin_orbit_ls_expect_from_coupling(int l, int j_2, int M_2);
```

```text
Direct verification of <L.S> using couple_states()'s CG coefficients,
evaluated via L.S = Lz Sz + (1/2) * (L+ S- + L- S+) acting term-by-term on
uncoupled product basis.
NOTE: This does not call closed-form spin_orbit_ls_expect() - it recomputes
<L.S> from scratch through coupling coefficients and l_plus_op/l_minus_op, as
independent cross-check that couple_states() and ladder operators agree with
analytic j(j+1)-l(l+1)-s(s+1) result.

Returns <L.S>/\hbar^2 (dimensionless), or NAN if (l,j_2,M_2) is invalid.
```
