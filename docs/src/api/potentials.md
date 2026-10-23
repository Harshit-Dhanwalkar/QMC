# `potentials`

> Generated from `physics/potentials.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Standard 1D potentials (infinite and finite well, harmonic, step, barrier, Morse, Coulomb, Yukawa) as V(x) functions.

**Category:** Single-particle quantum mechanics

```c
#include <qmc/physics/potentials.h>
```

**Start with:** [`V_harmonic`](#v_harmonic), [`V_finite_well`](#v_finite_well), [`V_barrier`](#v_barrier)

- Guide: [physics/potentials_1d.md](../physics/potentials_1d.md)
- Example: `examples/eg_05_uncertainty.c`
- Example: `examples/eg_06_hydrogen.c`
- Example: `examples/eg_07_central_potential.c`
- Test (every call with expected results): `tests/test_potentials.c`
- Runtime help: `qmc_help(stdout, "potentials");` or `qmc-help potentials`

## Types

### `potential_fn`

```text
Potential function pointer
```

```c
typedef double (*potential_fn)(double x, void *params);
```

## Functions (9)

[`V_infinite_well`](#v_infinite_well), [`V_finite_well`](#v_finite_well), [`V_harmonic`](#v_harmonic), [`V_step`](#v_step), [`V_barrier`](#v_barrier), [`V_coulomb`](#v_coulomb), [`V_yukawa`](#v_yukawa), [`V_morse`](#v_morse), [`potential_array`](#potential_array)

### `V_infinite_well`

```c
double V_infinite_well(double x, void *params);
```

```text
Infinite square well: V=0 for |x|<a, V=\infty elsewhere
```

### `V_finite_well`

```c
double V_finite_well(double x, void *params);
```

```text
Finite square well: V=-V0 for |x|<a, V=0 elsewhere
```

### `V_harmonic`

```c
double V_harmonic(double x, void *params);
```

```text
Harmonic oscillator: V = (1/2) * m * \omega^2 * x^2
```

### `V_step`

```c
double V_step(double x, void *params);
```

```text
Step potential: V=0 for x<0, V=V0 for x>=0
```

### `V_barrier`

```c
double V_barrier(double x, void *params);
```

```text
Rectangular barrier: V=V0 for 0<x<a, V=0 elsewhere
```

### `V_coulomb`

```c
double V_coulomb(double r, void *params);
```

```text
Coulomb potential: V = -k/r (for hydrogen)
```

### `V_yukawa`

```c
double V_yukawa(double r, void *params);
```

```text
Yukawa potential: V = -g * \exp(-\mu * r) / r (screening)
```

### `V_morse`

```c
double V_morse(double x, void *params);
```

```text
Morse potential: V = D*[1 - \exp(-a(x - x0))]^2 (molecules)
```

### `potential_array`

```c
void potential_array(double *x, int n, potential_fn V, void *params,
                     double *V_out);
```

```text
Generic: evaluate array of V values
```
