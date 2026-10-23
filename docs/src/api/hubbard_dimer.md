# `hubbard_dimer`

> Generated from `physics/hubbard_dimer.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Two-site Hubbard model: exact diagonalisation, eigenstates and time evolution.

**Category:** Condensed matter & many-body

```c
#include <qmc/physics/hubbard_dimer.h>
```

**Start with:** [`hd_solve`](#hd_solve), [`hd_n2_states`](#hd_n2_states), [`hd_basis_state`](#hd_basis_state), [`hd_eigenstate`](#hd_eigenstate), [`hd_evolve`](#hd_evolve)

- Example: `examples/eg_92_hubbard_dimer.c`
- Test (every call with expected results): `tests/test_hubbard_dimer.c`
- Runtime help: `qmc_help(stdout, "hubbard_dimer");` or `qmc-help hubbard_dimer`

## Overview

```text
Two-site Hubbard model (Hubbard dimer), solved by exact diagonalisation

  H = -t sum_{s} (c+_{1s} c_{2s} + h.c.) + U (n_{1up} n_{1dn} + n_{2up}
n_{2dn})

16-dimensional Fock space is built with Jordan-Wigner signs. Mode bits: 0 =
site 1 up, 1 = site 1 down, 2 = site 2 up, 3 = site 2 down; a basis state is
4-bit mask of occupied modes

At half filling (N = 2) spectrum is E = U/2 -+ sqrt(U^2/4 + 4 t^2) (two
singlets), U (antisymmetric charge-transfer singlet, which does not hop) and
0 (triplet, three-fold), so ground state is a singlet that crosses over from
a band-like state (U = 0) to a Heisenberg pair with J = 4 t^2 / U (U >> t)
A state of definite N and S_z has a diagonal one-site reduced density matrix,
whose entropy measures charge and spin entanglement between two sites

States are complex amplitude arrays re[16], im[16] indexed by mask
```

## Types

### `?`

```c
enum { HD_DIM = 16 };
```

### `hd_solution_t`

```c
typedef struct {
  double t, u;
  double e[HD_DIM];         /* energies, ascending */
  double v[HD_DIM][HD_DIM]; /* v[k][mask]: amplitude of basis state `mask` in
                               eigenstate k */
  int n[HD_DIM];            /* particle number of eigenstate k */
  int twosz[HD_DIM];        /* 2 S_z of eigenstate k */
} hd_solution_t;
```

### `hd_thermal_t`

```c
typedef struct {
  double energy, double_occ, ss,
      entropy; /* <H>, <n1up n1dn + n2up n2dn>, <S1.S2>, thermodynamic entropy S
                  / k_B */
} hd_thermal_t;
```

## Functions (12)

[`hd_solve`](#hd_solve), [`hd_n2_states`](#hd_n2_states), [`hd_basis_state`](#hd_basis_state), [`hd_eigenstate`](#hd_eigenstate), [`hd_evolve`](#hd_evolve), [`hd_expect_energy`](#hd_expect_energy), [`hd_expect_double_occ`](#hd_expect_double_occ), [`hd_expect_sz`](#hd_expect_sz), [`hd_expect_ss`](#hd_expect_ss), [`hd_site_probs`](#hd_site_probs), [`hd_site_entropy`](#hd_site_entropy), [`hd_thermal_n2`](#hd_thermal_n2)

### `hd_solve`

```c
int hd_solve(double t, double u, hd_solution_t *s);
```

```text
Diagonalises H; 0 on success, -1 for non-finite parameters
```

### `hd_n2_states`

```c
void hd_n2_states(const hd_solution_t *s, int idx[6]);
```

```text
Fills idx with six N = 2 eigenstates in ascending energy order (idx[0] is
ground state)
```

### `hd_basis_state`

```c
void hd_basis_state(int mask, double *re, double *im);
```

```text
Basis state |mask> as (re, im).
```

### `hd_eigenstate`

```c
void hd_eigenstate(const hd_solution_t *s, int k, double *re, double *im);
```

```text
Eigenstate k as (re, im).
```

### `hd_evolve`

```c
int hd_evolve(const hd_solution_t *s, const double *re0, const double *im0,
              double time, double *re, double *im);
```

```text
exp(-i H time) |psi0> (hbar = 1); -1 on non-finite time.
```

### `hd_expect_energy`

```c
double hd_expect_energy(const hd_solution_t *s, const double *re,
                        const double *im);
```

_No description in header yet._

### `hd_expect_double_occ`

```c
double hd_expect_double_occ(const double *re,
                            const double *im);
```

_No description in header yet._

### `hd_expect_sz`

```c
double hd_expect_sz(const double *re, const double *im,
                    int site);
```

```text
<n1up n1dn + n2up n2dn>
```

### `hd_expect_ss`

```c
double hd_expect_ss(const double *re, const double *im);
```

```text
<(n_up - n_dn) / 2> of site 0 or 1
```

### `hd_site_probs`

```c
void hd_site_probs(const double *re, const double *im, int site, double p[4]);
```

```text
Occupation probabilities of one site: p[0] empty, p[1] up, p[2] down, p[3]
double
```

### `hd_site_entropy`

```c
double hd_site_entropy(const double *re, const double *im, int site);
```

```text
Entropy -sum p ln p of those probabilities: entanglement entropy of a site
for a state of definite N and S_z
```

### `hd_thermal_n2`

```c
int hd_thermal_n2(const hd_solution_t *s, double temp, hd_thermal_t *out);
```

```text
Canonical (N = 2) thermal averages at temperature temp (k_B = 1); -1 unless
temp > 0
```
