# `quantum_walk`

> Generated from `physics/quantum_walk.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Discrete-time quantum walk on a line with a two-state coin.

**Category:** Condensed matter & many-body

```c
#include <qmc/physics/quantum_walk.h>
```

**Start with:** [`qwalk_init`](#qwalk_init), [`qwalk_step`](#qwalk_step), [`qwalk_norm`](#qwalk_norm), [`qwalk_probability`](#qwalk_probability), [`qwalk_mean`](#qwalk_mean)

- Guide: [physics/quantum_walk.md](../physics/quantum_walk.md)
- Example: `examples/eg_80_quantum_walk.c`
- Test (every call with expected results): `tests/test_quantum_walk.c`
- Runtime help: `qmc_help(stdout, "quantum_walk");` or `qmc-help quantum_walk`

## Overview

```text
Discrete-time quantum walk on a line

A walker with a two-state coin (up, down) at integer sites x = j - n/2
One step applies coin to every site, then shifts "up" one site to right and
"down" one site to left:

  [a']   [ \cos(th)  \sin(th)] [a]       a_new(x) = a'(x - 1)
  [b'] = [ \sin(th) -\cos(th)] [b]       b_new(x) = b'(x + 1)

th = \pi/4 is Hadamard walk. Map is unitary, so total probability is conserved
to round-off. Unlike classical random walk, whose width grows as sqrt(t),
quantum walk spreads ballistically: for a walker started with an unbiased
coin variance obeys

    <x^2> -> (1 - \sin(th)) t^2          (t >> 1)

and probability piles up near edges x = +/- t \cos(th) of light cone
("fronts"), with a lower density in middle. Classical walk, with same coin
flip probability 1/2 per step, has <x^2> = t

NOTE: Arrays have n sites and wrap around (periodic). Keep number of steps
below n/2 so probability from left edge has not reached right edge: walker
can move one site per step
```

## Macros and constants

```c
#define QWALK_COIN_UP 0   /* start in |up> */
```

```c
#define QWALK_COIN_DOWN 1 /* start in |down> */
```

```c
#define QWALK_COIN_SYM 2  /* start in (|up> + i |down>)/\sqrt(2), symmetric */
```

## Functions (8)

[`qwalk_init`](#qwalk_init), [`qwalk_step`](#qwalk_step), [`qwalk_norm`](#qwalk_norm), [`qwalk_probability`](#qwalk_probability), [`qwalk_mean`](#qwalk_mean), [`qwalk_variance`](#qwalk_variance), [`qwalk_classical`](#qwalk_classical), [`qwalk_asymptotic_variance`](#qwalk_asymptotic_variance)

### `qwalk_init`

```c
int qwalk_init(cvector_t *up, cvector_t *down, int coin);
```

```text
Put walker at centre site (n/2) with chosen coin state

Returns 0 on success, -1 on NULL, mismatched or n < 4 lengths, or bad `coin`
```

### `qwalk_step`

```c
int qwalk_step(cvector_t *up, cvector_t *down, double theta, int steps);
```

```text
Advance `steps` steps with coin angle \theta (radians).

Returns 0 on success, -1 on invalid input (NULL, mismatched length,
non-finite theta, steps < 1), -2 on allocation failure. On failure state is
untouched
```

### `qwalk_norm`

```c
double qwalk_norm(const cvector_t *up, const cvector_t *down);
```

```text
Total probability sum(|up|^2 + |down|^2); NaN on invalid input
```

### `qwalk_probability`

```c
int qwalk_probability(const cvector_t *up, const cvector_t *down, double *prob);
```

```text
P(x_j) = |up_j|^2 + |down_j|^2 into prob[0..n-1]. Returns 0 or -1
```

### `qwalk_mean`

```c
double qwalk_mean(const cvector_t *up, const cvector_t *down);
```

```text
Mean position and variance (in sites, centre site = 0); NaN on invalid
input or zero norm
```

### `qwalk_variance`

```c
double qwalk_variance(const cvector_t *up, const cvector_t *down);
```

_No description in the header yet._

### `qwalk_classical`

```c
int qwalk_classical(double *prob, int n, int t);
```

```text
Classical unbiased random walk after t steps from centre: P(x) =
C(t, (t + x)/2) / 2^t for x of same parity as t. Fills prob[0..n-1]
(index j is site j - n/2)

Returns 0, or -1 for NULL prob, n < 4, t < 0 or t >= n/2
```

### `qwalk_asymptotic_variance`

```c
double qwalk_asymptotic_variance(double theta);
```

```text
Limit of <x^2>/t^2 for quantum walk, [1 - \sin(\theta)] (\theta in [0, pi/2]);
NaN for non-finite input
```
