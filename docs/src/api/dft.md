# `dft`

> Generated from `physics/dft.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Kohn-Sham LDA density functional theory for closed-shell atoms and ions.

**Category:** Atomic structure & perturbation theory

```c
#include <qmc/physics/dft.h>
```

**Start with:** [`lda_exchange_energy_density`](#lda_exchange_energy_density), [`lda_exchange_potential`](#lda_exchange_potential), [`lda_correlation_energy_density_pz81`](#lda_correlation_energy_density_pz81), [`lda_correlation_potential_pz81`](#lda_correlation_potential_pz81), [`lda_xc_energy_density`](#lda_xc_energy_density)

- Example: `examples/eg_44_dft_atoms.c`
- Example: `examples/eg_55_molecular_dft.c`
- Example: `examples/eg_56_molecular_pbe_dft.c`
- Test (every call with expected results): `tests/test_dft.c`
- Runtime help: `qmc_help(stdout, "dft");` or `qmc-help dft`

## Overview

```text
NOTE: Kohn-Sham DFT, Local Density Approximation (LDA), for closed-shell
atoms/ions whose occupied subshells are all s-type.
The direct DFT analogue of hartree_fock.c's RHF. The only physical difference
from RHF is replaced exact exchange operator: RHF uses nonlocal Fock exchange
integral K; this uses local LDA exchange-correlation potential V_xc(n(r))
below, which is a function only of density at that point (not an integral
operator), making Kohn-Sham matrix purely diagonal in the potential.

Exchange-correlation functional: Slater/Dirac exchange (exact for uniform
electron gas) + Perdew-Zunger 1981 (PZ81) correlation (Perdew & Zunger, Phys.
Rev. B 23, 5048 (1981)), parametrized from Ceperley-Alder QMC data for the
unpolarized (paramagnetic, zeta=0) uniform electron gas.
```

## Types

### `dft_result_t`

```c
typedef struct {
  int n_orbitals;           /* number of doubly-occupied s-orbitals */
  int N;                    /* radial grid size */
  double Z;                 /* nuclear charge used */
  double *orbital_energies; /* size n_orbitals, converged KS eigenvalues */
  cvector_t **orbitals;     /* size n_orbitals; u_k(r) = r * R_k(r) in .re,
                             * normalized: integral u_k(r)^2 dr = 1 */
  double total_energy;      /* converged KS total energy (Hartree) */
  double E_hartree;         /* classical Hartree (electron-electron
                             * repulsion) energy component, informational */
  double E_xc;              /* exchange-correlation energy component,
                             * informational */
  int iterations;
  int converged;
} dft_result_t;
```

## Functions (16)

[`lda_exchange_energy_density`](#lda_exchange_energy_density), [`lda_exchange_potential`](#lda_exchange_potential), [`lda_correlation_energy_density_pz81`](#lda_correlation_energy_density_pz81), [`lda_correlation_potential_pz81`](#lda_correlation_potential_pz81), [`lda_xc_energy_density`](#lda_xc_energy_density), [`lda_xc_potential`](#lda_xc_potential), [`pw92_correlation_energy_density`](#pw92_correlation_energy_density), [`pw92_correlation_potential`](#pw92_correlation_potential), [`pbe_exchange_energy_density`](#pbe_exchange_energy_density), [`pbe_exchange_potential`](#pbe_exchange_potential), [`pbe_correlation_energy_density`](#pbe_correlation_energy_density), [`pbe_correlation_potential`](#pbe_correlation_potential), [`pbe_xc_energy_density`](#pbe_xc_energy_density), [`pbe_xc_potential`](#pbe_xc_potential), [`dft_lda_atom_s_orbitals`](#dft_lda_atom_s_orbitals), [`dft_result_free`](#dft_result_free)

### `lda_exchange_energy_density`

```c
double
lda_exchange_energy_density(double n);
```

```text
Slater/Dirac LDA exchange, n in electrons/bohr^3 (atomic units).
```

### `lda_exchange_potential`

```c
double lda_exchange_potential(double n);
```

```text
eps_x(n): energy per electron
```

### `lda_correlation_energy_density_pz81`

```c
double lda_correlation_energy_density_pz81(double n);
```

```text
PZ81 LDA correlation, unpolarized.
```

### `lda_correlation_potential_pz81`

```c
double lda_correlation_potential_pz81(double n);
```

_No description in header yet._

### `lda_xc_energy_density`

```c
double lda_xc_energy_density(double n);
```

```text
Sums of the above (what actually enters the Kohn-Sham potential/energy).
```

### `lda_xc_potential`

```c
double lda_xc_potential(double n);
```

_No description in header yet._

### `pw92_correlation_energy_density`

```c
double pw92_correlation_energy_density(double n);
```

```text
PBE (Perdew-Burke-Ernzerhof 1996) GGA exchange-correlation functional.
Reference: Perdew, Burke & Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996).

NOTE: PBE correlation's enhancement-factor H(rs,t) is analytically
parametrized against Perdew-Wang 1992 (PW92) LDA correlation (Perdew & Wang,
Phys. Rev. B 45, 13244 (1992)) -> NOT PZ81 above; so it needs its own
LDA-correlation base, pw92_correlation_*, distinct from
lda_correlation_*_pz81.

Unpolarized (\zeta=0, restricted/closed-shell) only
 \sigma = |\grad n|^2 (atomic units).
All energy densities are per-electron (eps_xc(n, \sigma)); multiply by n to
get energy density per volume.
```

### `pw92_correlation_potential`

```c
double pw92_correlation_potential(double n);
```

```text
eps_c(n), PW92, zeta=0
```

### `pbe_exchange_energy_density`

```c
double pbe_exchange_energy_density(double n, double sigma);
```

```text
PBE exchange eps_x(n,sigma) and its two functional derivatives:
 vrho = d(n * eps_x) / dn,
 vsigma = d(n*eps_x)/dsigma.
```

### `pbe_exchange_potential`

```c
void pbe_exchange_potential(double n, double sigma, double *vrho,
                            double *vsigma);
```

_No description in header yet._

### `pbe_correlation_energy_density`

```c
double pbe_correlation_energy_density(double n, double sigma);
```

```text
PBE correlation eps_c(n,sigma) = eps_c^PW92(n) + H(n, \sigma), and its two
functional derivatives.
```

### `pbe_correlation_potential`

```c
void pbe_correlation_potential(double n, double sigma, double *vrho,
                               double *vsigma);
```

_No description in header yet._

### `pbe_xc_energy_density`

```c
double pbe_xc_energy_density(double n, double sigma);
```

```text
Sums of the above (what actually enters the GGA Kohn-Sham potential/energy).
```

### `pbe_xc_potential`

```c
void pbe_xc_potential(double n, double sigma, double *vrho, double *vsigma);
```

_No description in header yet._

### `dft_lda_atom_s_orbitals`

```c
dft_result_t *dft_lda_atom_s_orbitals(const double *r, int N, double Z,
                                      int n_orbitals, double mix, double tol,
                                      int max_iter);
```

```text
Run closed-shell KS-LDA SCF for an s-orbitals-only atom/ion. Same argument
conventions as hartree_fock_atom_s_orbitals so two are directly comparable
side by side on same grid/Z/n_orbitals: r (uniform grid, r[0]>0 recommended),
N, Z, n_orbitals (electron count = 2 *  n_orbitals, all in s subshells: He=1,
Be=2, ...), mix (linear density-mixing damping in (0,1]), tol (convergence
threshold on successive eigenvalue change), max_iter.

Returns NULL on invalid input / allocation failure.
```

### `dft_result_free`

```c
void dft_result_free(dft_result_t *res);
```

_No description in header yet._
