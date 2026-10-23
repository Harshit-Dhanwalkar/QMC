# `identical`

> Generated from `physics/identical.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Identical particles: Slater determinants for fermions and permanents for bosons.

**Category:** Atomic structure & perturbation theory

```c
#include <qmc/physics/identical.h>
```

**Start with:** [`slater_matrix`](#slater_matrix), [`slater_determinant_value`](#slater_determinant_value), [`bosonic_permanent_value`](#bosonic_permanent_value)

- Guide: [physics/identical.md](../physics/identical.md)
- Example: `examples/eg_09_identical_particles.c`
- Test (every call with expected results): `tests/test_identical.c`
- Runtime help: `qmc_help(stdout, "identical");` or `qmc-help identical`

## Overview

```text
Identical particles: build the N x N matrix needed to evaluate an
N-particle Slater determinant (fermions) or permanent (bosons) at a
specific configuration, and compute that determinant/permanent.

orbitals[i]: tabulated single-particle orbital i, sampled on a common
             position grid of length M.
indices[j] : grid index of particle j's position (0 <= indices[j] < M).
Two particles at same grid index (indices[j1]==indices[j2]) or two identical
orbitals both make Slater determinant vanish exactly (Pauli exclusion)
```

## Functions (3)

[`slater_matrix`](#slater_matrix), [`slater_determinant_value`](#slater_determinant_value), [`bosonic_permanent_value`](#bosonic_permanent_value)

### `slater_matrix`

```c
cmatrix_t *slater_matrix(cvector_t **orbitals, int N, const int *indices);
```

```text
M_ij = orbitals[i]->data[indices[j]], size N x N
```

### `slater_determinant_value`

```c
complex_t slater_determinant_value(cvector_t **orbitals, int N,
                                   const int *indices);
```

```text
Fermions: \psi(x_1,...,x_N) = (1 / \sqrt(N!)) * det[\phi_i(x_j)]
Exactly 0 if two particles share a position or two orbitals coincide.
```

### `bosonic_permanent_value`

```c
complex_t bosonic_permanent_value(cvector_t **orbitals, int N,
                                  const int *indices);
```

```text
Bosons: \psi(x_1,...,x_N) = (1 / \sqrt(N! * prod_k n_k!)) * perm[\phi_i(x_j)]
is invariant under exchanging any two particles. n_k is the occupation number
of each distinct orbital among orbitals[0..n-1]; this reduces to simpler
1 / sqrt(N!) exactly when every orbital is distinct, but prod_k n_k! factor
matters whenever multiple bosons share an orbital
```
