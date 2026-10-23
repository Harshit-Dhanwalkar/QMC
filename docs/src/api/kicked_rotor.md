# `kicked_rotor`

> Generated from `physics/kicked_rotor.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Quantum kicked rotor and the classical standard map.

**Category:** Two-level, open systems & quantum optics

```c
#include <qmc/physics/kicked_rotor.h>
```

**Start with:** [`krotor_init`](#krotor_init), [`krotor_step`](#krotor_step), [`krotor_momentum_distribution`](#krotor_momentum_distribution), [`krotor_momentum2`](#krotor_momentum2), [`krotor_norm`](#krotor_norm)

- Example: `examples/eg_82_kicked_rotor.c`
- Test (every call with expected results): `tests/test_kicked_rotor.c`
- Runtime help: `qmc_help(stdout, "kicked_rotor");` or `qmc-help kicked_rotor`

## Overview

```text
Quantum kicked rotor and classical standard map

  H = p^2/2 + K \cos(\theta) \sum_n \delta(t - n),   \theta in [0, 2\pi)

Classical: the Chirikov standard map
    p'     = p + K \sin(\theta)
    \theta' = \theta + p'           (mod 2 \pi)

is regular for small K and globally chaotic for K > 1 (K_c ~ 0.9716)
NOTE: In chaotic regime momentum diffuses, <p^2> ~ D n with D ~ K^2/2
(quasilinear value, which holds up to Bessel-function corrections)

Quantum: p = \hbar m with integer m, and one period of Floquet map is
    psi -> \exp(-i \hbar m^2 / 2) . \exp(-i (K / \hbar) cos(\theta)) \psi

applied as a kick in angle space and free rotation in momentum space (FFT
between them). Here \hbar is dimensionless effective Planck constant

Quantum interference stops diffusion after break time: <m^2> saturates
(dynamical localisation, momentum-space cousin of Anderson localisation) and
momentum distribution has exponential tails

Exact special case, quantum resonance \hbar = 4\pi: free factor is 1, so
after n kicks \psi = \exp(-i n (K/\hbar) \cos(\theta)) and momentum
distribution is |J_m(n K / \hbar)|^2, giving <m^2> = (n K/hbar)^2 / 2 for
every n (ballistic, not diffusive)

NOTE: state lives on N = 2^k angle points \theta_j = 2 \pi j / N and momentum
index m runs over -N/2 .. N/2 - 1; keep distribution well inside that range
```

## Functions (8)

[`krotor_init`](#krotor_init), [`krotor_step`](#krotor_step), [`krotor_momentum_distribution`](#krotor_momentum_distribution), [`krotor_momentum2`](#krotor_momentum2), [`krotor_norm`](#krotor_norm), [`krotor_classical_init`](#krotor_classical_init), [`krotor_classical_step`](#krotor_classical_step), [`krotor_classical_p2`](#krotor_classical_p2)

### `krotor_init`

```c
int krotor_init(cvector_t *psi);
```

```text
Start in momentum eigenstate m = 0 (uniform in angle)

Returns 0, or -1 for NULL or a length that is not a power of two >= 4
```

### `krotor_step`

```c
int krotor_step(cvector_t *psi, double k, double hbar, int steps);
```

```text
Apply `steps` kicks + free rotations with kick strength k and effective
Planck constant \hbar

Returns 0, or -1 for invalid input (NULL, bad length, non-finite k, hbar <= 0
or non-finite, steps < 1), -2 on allocation failure Norm is conserved to
round-off
```

### `krotor_momentum_distribution`

```c
int krotor_momentum_distribution(const cvector_t *psi, double *prob);
```

```text
Momentum probabilities |c_m|^2 into prob[0..N-1]; index j is m = j - N/2.
Returns 0, or -1 on invalid input, -2 on allocation failure
```

### `krotor_momentum2`

```c
double krotor_momentum2(const cvector_t *psi);
```

```text
<m^2> and norm; NaN on invalid input
```

### `krotor_norm`

```c
double krotor_norm(const cvector_t *psi);
```

_No description in the header yet._

### `krotor_classical_init`

```c
int krotor_classical_init(double *theta, double *p, int n, double hbar,
                          unsigned long long seed);
```

```text
Classical ensemble of n points for standard map. Initial angles are uniform
in [0, 2 pi) and momenta uniform in [-hbar/2, hbar/2) (width of quantum m = 0
state), from a deterministic generator seeded by `seed`

Returns 0, or -1 for NULL arrays or n < 1
```

### `krotor_classical_step`

```c
int krotor_classical_step(double *theta, double *p, int n, double k, int steps);
```

```text
Advance ensemble `steps` times (\theta wrapped into [0, 2 \pi); p is not
wrapped)

Returns 0, or -1 for NULL, n < 1, non-finite k or steps < 1
```

### `krotor_classical_p2`

```c
double krotor_classical_p2(const double *p, int n);
```

```text
Mean of p^2 over ensemble; NaN for NULL or n < 1
```
