# `uncertainty`

> Generated from `physics/uncertainty.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Position and momentum uncertainties and energy expectation values for wavefunctions.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/uncertainty.h>
```

**Start with:** [`compute_uncertainties`](#compute_uncertainties), [`compute_energy_expectation`](#compute_energy_expectation)

- Guide: [physics/uncertainty.md](../physics/uncertainty.md)
- Example: `examples/eg_05_uncertainty.c`
- Test (every call with expected results): `tests/test_uncertainty.c`
- Runtime help: `qmc_help(stdout, "uncertainty");` or `qmc-help uncertainty`

## Types

### `uncertainty_t`

```text
Container for uncertainty results
```

```c
typedef struct {
  double mean_x;
  double mean_p;
  double delta_x;
  double delta_p;
  double product; /* \Delta x * \Delta p */
} uncertainty_t;
```

## Functions (2)

[`compute_uncertainties`](#compute_uncertainties), [`compute_energy_expectation`](#compute_energy_expectation)

### `compute_uncertainties`

```c
uncertainty_t compute_uncertainties(const wavefunction_t *wf);
```

```text
Compute uncertainties for wavefunction
```

### `compute_energy_expectation`

```c
double compute_energy_expectation(const wavefunction_t *wf, potential_fn V,
                                  void *params, double mass);
```

```text
Energy expectation <H> for given potential and particle mass
```
