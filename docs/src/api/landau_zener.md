# `landau_zener`

> Generated from `physics/landau_zener.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Landau-Zener sweeps of a two-level system on the Bloch sphere.

**Category:** Two-level, open systems & quantum optics

```c
#include <qmc/physics/landau_zener.h>
```

**Start with:** [`lz_probability`](#lz_probability), [`lz_upper_population`](#lz_upper_population), [`lz_lower_state`](#lz_lower_state), [`lz_advance`](#lz_advance), [`lz_sweep`](#lz_sweep)

- Guide: [physics/landau_zener.md](../physics/landau_zener.md)
- Example: `examples/eg_81_landau_zener.c`
- Test (every call with expected results): `tests/test_landau_zener.c`
- Runtime help: `qmc_help(stdout, "landau_zener");` or `qmc-help landau_zener`

## Overview

```text
Landau-Zener sweeps of a two-level system on Bloch sphere (\hbar = 1)

  H(t)  = (1/2) [ \Omega \sigma_x + \delta(t) \sigma_z ]
  dv/dt = b x v
  b     = (\Omega, 0, \delta)

(Convention: Bloch vector v, basis |g> = north pole)
Detuning is swept linearly through resonance at rate `rate` (d \delta/dt) so
2 diabatic levels cross at \delta = 0 and adiabatic levels (+/-
\sqrt(\Omega^2 + \delta^2)/2) avoid each other with a minimum gap \Omega

Landau-Zener formula: a sweep from delta = -infinity to +infinity that starts
in lower adiabatic state ends in upper one (a diabatic jump, i.e. state keeps
its diabatic character) with probability

    P_LZ = \exp(-\pi \Omega^2 / (2 |rate|))

Sweeping out and back (two passes) makes two paths interfere, giving
Stueckelberg oscillations bounded by 4 P_LZ (1 - P_LZ)
```

## Functions (5)

[`lz_probability`](#lz_probability), [`lz_upper_population`](#lz_upper_population), [`lz_lower_state`](#lz_lower_state), [`lz_advance`](#lz_advance), [`lz_sweep`](#lz_sweep)

### `lz_probability`

```c
double lz_probability(double omega, double rate);
```

```text
Closed-form probability exp(-pi omega^2 / (2 |rate|)); NaN for invalid input
(non-finite omega, rate == 0 or non-finite)
```

### `lz_upper_population`

```c
double lz_upper_population(const double v[3], double omega, double delta);
```

```text
Population of upper adiabatic level, (1 + v . b/|b|)/2 for field
b = (\omega, 0, \delta); NaN for non-finite input or b = 0
```

### `lz_lower_state`

```c
int lz_lower_state(double v[3], double omega, double delta);
```

```text
Put v in lower adiabatic eigenstate for given field,
v = -b/|b|

Returns 0, or -1 for NULL v, non-finite input or b = 0
```

### `lz_advance`

```c
int lz_advance(double v[3], double omega, double delta0, double rate,
               double duration, double dt);
```

```text
Integrate dv/dt = b(t) x v with RK4 for `duration` while detuning moves
linearly, \delta(t) = \delta_0 + rate t (rate may be negative or zero)
Steps are at most dt and at most 0.03 / sqrt(omega^2 + dmax^2), dmax largest
|\delta| reached, so |v| is conserved to ~1e-9 even though RK4 is not exactly
norm preserving. v is updated in place

Returns 0 on success, -1 on invalid input (NULL v, non-finite or |v| > 1 +
1e-9, \omega < 0 or non-finite, non-finite delta0 / rate, duration < 0 or
non-finite, dt <= 0, more than 100 million steps). On failure v is untouched
```

### `lz_sweep`

```c
int lz_sweep(double v[3], double omega, double rate, double amp, double dt,
             int passes);
```

```text
Integrate dv/dt = b(t) x v with RK4 while detuning zigzags between -amp and
+amp at speed |rate|: starting at delta = -amp, ramp up to +amp (pass 1),
then down to -amp (pass 2), and so on for `passes` passes. Each pass takes 2
amp/|rate| and is integrated with lz_advance. v is updated in place; for
standard Landau-Zener experiment start with lz_lower_state(v, omega, -amp)

Returns 0 on success, -1 on invalid input (NULL v, non-finite or |v| > 1 +
1e-9, omega < 0 or non-finite, rate == 0 or non-finite, amp <= 0 or
non-finite, dt <= 0, passes < 1, more than 100 million steps). On failure v
is untouched
```
