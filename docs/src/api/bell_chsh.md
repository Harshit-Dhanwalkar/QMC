# `bell_chsh`

> Generated from `physics/bell_chsh.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Bell and CHSH correlations: quantum values and local-hidden-variable and deterministic bounds.

**Category:** Quantum computing & information

```c
#include <qmc/physics/bell_chsh.h>
```

**Start with:** [`bell_correlation`](#bell_correlation), [`bell_chsh`](#bell_chsh), [`bell_local_correlation`](#bell_local_correlation), [`bell_local_chsh`](#bell_local_chsh), [`bell_deterministic_bound`](#bell_deterministic_bound)

- Example: `examples/eg_88_bell_chsh.c`
- Test (every call with expected results): `tests/test_bell_chsh.c`
- Runtime help: `qmc_help(stdout, "bell_chsh");` or `qmc-help bell_chsh`

## Overview

```text
CHSH Bell test on a two-qubit Werner state

\rho = V |\psi-><\psi-| + (1 - V) I / 4, singlet mixed with white noise
(visibility V). Alice and Bob measure spin along directions in x-z plane
at angles a and b, A(a) = \cos(a) \sigma_z + \sin(a) \sigma_x. Their
correlation is E(a, b) = Tr[\rho A(a) (x) A(b)] = -V \cos(a - b), computed
here from density matrix, not from closed form

CHSH quantity S = E(a,b) - E(a,b') + E(a',b) + E(a',b'). Any local hidden
variable theory has |S| <= 2; quantum mechanics reaches 2 \sqrt{2} V
(Tsirelson), so Bell inequality is violated for V > 1 / \sqrt{2}
A concrete local model (a hidden angle, each side outputs sign of projection)
gives sawtooth correlation bell_local_correlation
```

## Functions (7)

[`bell_correlation`](#bell_correlation), [`bell_chsh`](#bell_chsh), [`bell_local_correlation`](#bell_local_correlation), [`bell_local_chsh`](#bell_local_chsh), [`bell_deterministic_bound`](#bell_deterministic_bound), [`bell_max_chsh`](#bell_max_chsh), [`bell_sample`](#bell_sample)

### `bell_correlation`

```c
double bell_correlation(double visibility, double a, double b);
```

```text
E(a, b) from density matrix; NaN unless 0 <= V <= 1 and angles are finite
```

### `bell_chsh`

```c
double bell_chsh(double visibility, double a, double a2, double b, double b2);
```

```text
S = E(a,b) - E(a,b2) + E(a2,b) + E(a2,b2); NaN on bad input
```

### `bell_local_correlation`

```c
double bell_local_correlation(double a, double b);
```

```text
Correlation of hidden-angle local model: -(1 - 2 d / \pi), d = |a - b| folded
into [0, \pi]
```

### `bell_local_chsh`

```c
double bell_local_chsh(double a, double a2, double b, double b2);
```

```text
CHSH value of that local model for same settings
```

### `bell_deterministic_bound`

```c
double bell_deterministic_bound(void);
```

```text
Largest |S| over 16 deterministic local strategies (always 2)
```

### `bell_max_chsh`

```c
double bell_max_chsh(double visibility, int n_grid);
```

```text
Largest |S| found by scanning a 4-angle grid with n_grid steps per angle over
[0, 2 \pi); NaN on bad input
```

### `bell_sample`

```c
double bell_sample(double visibility, double a, double a2, double b, double b2,
                   int n_pairs, unsigned seed, double outcome_sums[4]);
```

```text
Simulate n_pairs measurements for each of four setting pairs (a,b), (a,b2),
(a2,b), (a2,b2) with a deterministic generator (seed)
outcome_sums[i] gets sum of products A*B (so E_est = outcome_sums[i] /
n_pairs)

Returns estimate of S, or NaN on bad input
```
