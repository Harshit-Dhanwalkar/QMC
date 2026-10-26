# `angular`

> Generated from `physics/angular.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Angular momentum and spin: Pauli matrices, L+/L-/L^2/Lz, Clebsch-Gordan and Gaunt coefficients.

**Category:** Atomic structure & perturbation theory

```c
#include <qmc/physics/angular.h>
```

**Start with:** [`lz_matrix`](#lz_matrix), [`l2_matrix`](#l2_matrix), [`clebsch_gordan`](#clebsch_gordan)

- Guide: [physics/angular_momentum.md](../physics/angular_momentum.md)
- Example: `examples/eg_15_angular_coupling.c`
- Test (every call with expected results): `tests/test_angular_coupling.c`
- Runtime help: `qmc_help(stdout, "angular");` or `qmc-help angular`

## Data

```c
extern complex_t sigma_x[4];
```

```c
extern complex_t sigma_y[4];
```

```c
extern complex_t sigma_z[4];
```

## Functions (13)

[`l_plus_op`](#l_plus_op), [`l_minus_op`](#l_minus_op), [`l2_matrix`](#l2_matrix), [`lz_matrix`](#lz_matrix), [`lx_matrix`](#lx_matrix), [`ly_matrix`](#ly_matrix), [`clebsch_gordan`](#clebsch_gordan), [`couple_states`](#couple_states), [`couple_allowed_J`](#couple_allowed_j), [`gaunt_coefficient`](#gaunt_coefficient), [`spin_sigma_x`](#spin_sigma_x), [`spin_sigma_y`](#spin_sigma_y), [`spin_sigma_z`](#spin_sigma_z)

### `l_plus_op`

```c
complex_t l_plus_op(int orbital_l, int proj_m,
                    int proj_m_prime);
```

```text
Ladder operators for angular momentum: L_{\pm} = L_x \pm i L_y
```

### `l_minus_op`

```c
complex_t l_minus_op(int orbital_l, int proj_m,
                     int proj_m_prime);
```

```text
matrix element <l,m+1|L_+|l,m>
```

### `l2_matrix`

```c
cmatrix_t *l2_matrix(int l);
```

```text
Matrix representation of L^2, L_z for a given l (dimension 2l+1)
```

### `lz_matrix`

```c
cmatrix_t *lz_matrix(int l);
```

_No description in header yet._

### `lx_matrix`

```c
cmatrix_t *lx_matrix(int l);
```

_No description in header yet._

### `ly_matrix`

```c
cmatrix_t *ly_matrix(int l);
```

_No description in header yet._

### `clebsch_gordan`

```c
double clebsch_gordan(int j1_2, int m1_2, int j2_2, int m2_2, int J_2, int M_2);
```

```text
Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M>, general j1/j2/J via
Racah's formula

  spin-1/2, m=+1/2    -> j=1, m=+1
  orbital   l=1, m=0  -> j=2, m=0
  total     J=3/2     -> J=3

Returns 0 if inputs violate triangle inequality, m1+m2 != M, |m| > j for any
of three, or required integer/half-integer parity (j1+j2+J must be an
integer)
Phase convention: Condon-Shortley
```

### `couple_states`

```c
cvector_t *couple_states(int j1_2, int j2_2, int J_2, int M_2);
```

```text
Angular momentum coupling: |j1,j2; J,M> in uncoupled product basis {|j1,m1>
(x) |j2,m2>}
```

### `couple_allowed_J`

```c
int couple_allowed_J(int j1_2, int j2_2, int *J2_out);
```

```text
Enumerate allowed total-J values (doubled) for given j1_2, j2_2, i.e.
J_2 = |j1_2-j2_2|, |j1_2-j2_2|+2, ..., j1_2+j2_2
```

### `gaunt_coefficient`

```c
double gaunt_coefficient(int l, int m, int k, int lp, int mp);
```

```text
Gaunt coefficient c^k(l,m; l',m') = \sqrt(4 * \pi / (2k+1)) * \int Y_lm^*
Y_kq Y_l'm' d\Omega, q = m - m'. Built from clebsch_gordan() above via
  c^k(l,m;l',m') = (-1)^(m-m') * sqrt((2l+1)/(2l'+1))
                   * <l0;k0|l'0> * <l,-m;k,q|l',-m'>
l, m, k, l', m' are integers, since orbital angular momenta are always
integers

Returns 0 if triangle inequality |l-l'| <= k <= l+l' is violated, if l+k+l'
is odd or if |m|>l or |m'|>l'
```

### `spin_sigma_x`

```c
void spin_sigma_x(cvector_t *spinor);
```

```text
Spin-1/2 operations: apply \sigma_x, \sigma_y, \sigma_z to 2-component spinor
```

### `spin_sigma_y`

```c
void spin_sigma_y(cvector_t *spinor);
```

_No description in header yet._

### `spin_sigma_z`

```c
void spin_sigma_z(cvector_t *spinor);
```

_No description in header yet._
