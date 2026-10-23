# `grover_multi`

> Generated from `physics/grover_multi.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Grover search with several marked items: state evolution and success probability.

**Category:** Quantum computing & information

```c
#include <qmc/physics/grover_multi.h>
```

**Start with:** [`grover_multi_init`](#grover_multi_init), [`grover_multi_step`](#grover_multi_step), [`grover_multi_success`](#grover_multi_success), [`grover_multi_theta`](#grover_multi_theta), [`grover_multi_success_exact`](#grover_multi_success_exact)

- Example: `examples/eg_87_grover_multi.c`
- Test (every call with expected results): `tests/test_grover_multi.c`
- Runtime help: `qmc_help(stdout, "grover_multi");` or `qmc-help grover_multi`

## Overview

```text
Grover search with several marked items, on full state vector

N = 2^n items, M of them marked. One Grover iteration is oracle (flip sign of
every marked amplitude) followed by diffusion operator 2|s><s| - 1 (inversion
about mean). State stays in plane spanned by uniform marked and unmarked
superpositions and rotates by 2 \theta per iteration, \sin(\theta) = \sqrt(M
/ N), so success probability after k iterations is \sin^2((2k + 1) \theta)
and best k is near \pi / (4 \theta) - 1/2

State is real: amplitudes are stored as doubles
```

## Macros and constants

```c
#define GROVER_MAX_QUBITS 12
```

## Functions (6)

[`grover_multi_init`](#grover_multi_init), [`grover_multi_step`](#grover_multi_step), [`grover_multi_success`](#grover_multi_success), [`grover_multi_theta`](#grover_multi_theta), [`grover_multi_success_exact`](#grover_multi_success_exact), [`grover_multi_optimal_iterations`](#grover_multi_optimal_iterations)

### `grover_multi_init`

```c
int grover_multi_init(double *amp, int n_qubits);
```

```text
amp[0 .. 2^n - 1] = uniform superposition 1 / \sqrt(N). 0 on success, -1 on
bad input
```

### `grover_multi_step`

```c
int grover_multi_step(double *amp, int n_qubits, const unsigned char *marked);
```

```text
One iteration; marked[i] != 0 flags a marked item. 0 on success, -1 on bad
input
```

### `grover_multi_success`

```c
double grover_multi_success(const double *amp, int n_qubits,
                            const unsigned char *marked);
```

```text
Total probability on marked items; NaN on bad input
```

### `grover_multi_theta`

```c
double grover_multi_theta(int n_qubits, int n_marked);
```

```text
Rotation angle \theta = \asin(\sqrt(M / N)); NaN unless 1 <= M <= N (M = N
gives \pi / 2)
```

### `grover_multi_success_exact`

```c
double grover_multi_success_exact(int n_qubits, int n_marked, int k);
```

```text
Closed-form success probability \sin^2((2k + 1) \theta)
```

### `grover_multi_optimal_iterations`

```c
int grover_multi_optimal_iterations(int n_qubits, int n_marked);
```

```text
Iteration count with highest success probability, round(\pi / (4 \theta) -
1/2); -1 on bad input
```
