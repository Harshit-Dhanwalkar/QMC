# `variational`

> Generated from `physics/variational.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Variational method: generic golden-section minimiser and variational energy minimisation.

**Category:** Atomic structure & perturbation theory

```c
#include <qmc/physics/variational.h>
```

**Start with:** [`golden_section_minimize`](#golden_section_minimize), [`variational_energy`](#variational_energy), [`variational_minimize`](#variational_minimize)

- Guide: [physics/variational.md](../physics/variational.md)
- Example: `examples/eg_73_floquet.c`
- Test (every call with expected results): `tests/test_floquet.c`
- Runtime help: `qmc_help(stdout, "variational");` or `qmc-help variational`

## Types

### `trial_func`

```c
typedef struct {
  void (*trial_func)(double alpha, wavefunction_t *wf);
  wavefunction_t *wf;
  potential_fn V;
  void *params;
  double mass;
} variational_closure_t;
```

## Functions (3)

[`golden_section_minimize`](#golden_section_minimize), [`variational_energy`](#variational_energy), [`variational_minimize`](#variational_minimize)

### `golden_section_minimize`

```c
double golden_section_minimize(double a, double b,
                               double (*f)(double x, void *params),
                               void *params, double tol);
```

```text
Generic 1D golden-section minimizer of arbitrary scalar function
f(x, params) over [a,b].
```

### `variational_energy`

```c
double variational_energy(const wavefunction_t *wf, potential_fn V,
                          void *params, double mass);
```

```text
Compute expectation value of Hamiltonian for trial wavefunction
```

### `variational_minimize`

```c
double variational_minimize(double alpha_min, double alpha_max,
                            void (*trial_func)(double alpha,
                                               wavefunction_t *wf),
                            wavefunction_t *wf, potential_fn V, void *params,
                            double mass, double tol);
```

```text
Minimize energy with respect to single parameter alpha
```
