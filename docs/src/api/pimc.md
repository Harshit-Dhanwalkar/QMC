# `pimc`

> Generated from `physics/pimc.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Path integral Monte Carlo for two-electron atoms at finite temperature.

**Category:** Quantum Monte Carlo

```c
#include <qmc/physics/pimc.h>
```

**Start with:** [`kelbg_potential`](#kelbg_potential), [`kelbg_energy_correction`](#kelbg_energy_correction), [`pimc_walker_alloc`](#pimc_walker_alloc), [`pimc_walker_init`](#pimc_walker_init), [`pimc_bisection_move`](#pimc_bisection_move)

- Guide: [physics/pimc.md](../physics/pimc.md)
- Example: `examples/eg_30_pimc_helium.c`
- Test (every call with expected results): `tests/test_pimc.c`
- Runtime help: `qmc_help(stdout, "pimc");` or `qmc-help pimc`

## Overview

```text
Path Integral Monte Carlo (PIMC) for two-electron atom/ion (nuclear charge Z,
fixed at origin), finite-temperature imaginary-time path integral method

2 electrons are treated as distinguishable ("boltzmannons")

Coulomb singularity: bare Coulomb potential is singular at r=0, and naive
discretized (primitive or higher-order) action applied directly to it
produces a "path collapse" catastrophe when a bead's random walk brings two
charges close together. This module instead uses Kelbg pair potential (
Reference: Kelbg 1961; see e.g. Filinov et al., "Dynamical Properties and
Plasmon Dispersion of Weakly Degenerate Correlated One-Component Plasma",
arXiv:physics/0012027, eq. 10), which is finite at r=0 and reduces to bare
Coulomb potential in \tau->0 (high-temperature) limit

Sampling: bisection moves (Ceperley 1995), A full-ring bisection resamples an
entire electron's P-bead path in one proposal, drawn exactly from
free-particle (kinetic) bridge distribution, so Metropolis acceptance depends
only on potential-action difference

All energies in Hartree atomic units (hbar = m_e = 1)
```

## Types

### `pimc_walker_t`

```c
typedef struct {
  int P;           /* number of imaginary-time beads (ring polymer length),
                      must be a power of 2 */
  double (*r1)[3]; /* electron 1 beads, r1[0..P-1] */
  double (*r2)[3]; /* electron 2 beads, r2[0..P-1] */
} pimc_walker_t;
```

### `pimc_result_t`

```c
typedef struct {
  double energy; /* block-averaged thermodynamic energy estimator, Hartree */
  double error;  /* standard error on energy */
  double energy_virial;   /* block-averaged virial energy estimator, Hartree */
  double error_virial;    /* standard error on energy_virial */
  double acceptance_rate; /* fraction of bisection moves accepted */
  int n_blocks;           /* number of statistics blocks completed */
} pimc_result_t;
```

## Functions (10)

[`kelbg_potential`](#kelbg_potential), [`kelbg_energy_correction`](#kelbg_energy_correction), [`pimc_walker_alloc`](#pimc_walker_alloc), [`pimc_walker_free`](#pimc_walker_free), [`pimc_walker_init`](#pimc_walker_init), [`pimc_bisection_move`](#pimc_bisection_move), [`pimc_energy_estimator`](#pimc_energy_estimator), [`pimc_virial_estimator`](#pimc_virial_estimator), [`pimc_run`](#pimc_run), [`pimc_run_parallel`](#pimc_run_parallel)

### `kelbg_potential`

```c
double kelbg_potential(double r, double lambda, double q);
```

```text
Kelbg-regularized pair potential between two charges q_a, q_b (in atomic
units, q_a*q_b so sign of product gives attraction/repulsion) separated by
distance r, at imaginary-time-slice thermal wavelength:
   \lambda = \sqrt(\tau / \mu)
Where
  \mu = pair reduced mass

  V_Kelbg(r) = (q/r) * (1 - \exp(-r^2 / \lambda^2)) + q * (\sqrt(\pi) /
  \lambda) * erfc(r / \lambda),
  q = q_a*q_b

Finite at r=0 (V_Kelbg(0) = q * \sqrt(\pi) / \lambda) reduces to bare Coulomb
q/r as \lambda -> 0 (\tau -> 0).
```

### `kelbg_energy_correction`

```c
double kelbg_energy_correction(double r, double lambda, double q);
```

```text
V_Kelbg(r, \tau) + \tau * dV_Kelbg/ d\tau, potential-side contribution to
thermodynamic energy estimator when using a tau-dependent effective
potential. Closed form:
  V_Kelbg(r, \tau) + \tau * dV / d\tau = (q/r) * (1 - \exp(-r^2 / \lambda^2))
+ q * (\sqrt(\pi) / (2 * \lambda)) * erfc(r / \lambda)

i.e. identical to V_Kelbg except erfc term's coefficient is halved.
```

### `pimc_walker_alloc`

```c
pimc_walker_t *pimc_walker_alloc(int P);
```

_No description in the header yet._

### `pimc_walker_free`

```c
void pimc_walker_free(pimc_walker_t *w);
```

_No description in the header yet._

### `pimc_walker_init`

```c
void pimc_walker_init(pimc_walker_t *w, rng_state_t *rng, double Z);
```

```text
Initialize both electrons' P beads to small random cluster around (+/-0.5, 0,
0)-ish
```

### `pimc_bisection_move`

```c
int pimc_bisection_move(pimc_walker_t *w, int which, double Z, double tau,
                        int level, rng_state_t *rng);
```

```text
Bisection move for electron `which` (0 or 1): resamples an interior segment
of length 2^level beads from exact free-particle bridge distribution, then
accepts/rejects via Metropolis on potential-action difference over just that
segment

Returns 1 if move was accepted, 0 if rejected
```

### `pimc_energy_estimator`

```c
double pimc_energy_estimator(const pimc_walker_t *w, double Z, double tau);
```

```text
One instantaneous energy estimator sample (Hartree) at walker's current
configuration
NOTE: Callers should average this over many post-equilibration configurations
```

### `pimc_virial_estimator`

```c
double pimc_virial_estimator(const pimc_walker_t *w, double Z, double tau);
```

```text
One instantaneous virial energy estimator sample (Hartree) at walker's
current configuration

NOTE: An alternative to pimc_energy_estimator with same expectation value but
substantially lower variance (thermodynamic estimator's kinetic term is a
difference of two large, nearly-cancelling quantities that grows noisier as P
increases, virial estimator sidesteps this entirely by eliminating explicit
kinetic-energy difference term)

Ring-polymer partition function Z(\beta) is invariant under relabeling
beads' displacements from each particle's centroid R_c by any scale
factor \lambda, differentiating both sides at \lambda=1 gives an exact
relation between kinetic-energy term and potential-gradient term:

  E_virial = dN/(2 * \beta)
             + (1/P) \sum_i [ V(R_i, \tau) + \tau * dV / d\tau(R_i, \tau)
                      + (1/2) \sum_p (r_i^p - r_c^p).grad_p V(R_i, \tau) ]

(\tau-fixed) spatial gradient:
  dV_Kelbg/dr = -q * (1 - \exp(-r^2 / \lambda^2)) / r^2
```

### `pimc_run`

```c
pimc_result_t pimc_run(double Z, int P, double tau, int level,
                       int n_equilibration, int n_blocks, int block_size,
                       uint64_t seed);
```

```text
Run PIMC for two-electron atom/ion of nuclear charge Z: P beads,
imaginary-time step \tau (so \beta = P * \tau; larger \beta = lower
temperature, needed for estimator to approach ground-state energy).
Equilibrates for n_equilibration sweeps (discarded; one sweep = P/2^level
bisection moves per electron, covering ring once), then runs n_blocks
statistics blocks of block_size sweeps each, measuring energy estimator once
per sweep and block-averaging for reported energy an error

NOTE: P must be a power of 2. level sets bisection segment length (2^level
beads); must satisfy 1 <= level <= log2(P). Smaller segments keep
acceptance high at large P (a full-ring move, level=log2(P), has
acceptance collapsing toward 0 as P grows, since whole path is resampled in
one proposal); level in 3-6 range is a reasonable starting point for most P
```

### `pimc_run_parallel`

```c
pimc_result_t pimc_run_parallel(int n_replicas, double Z, int P, double tau,
                                int level, int n_equilibration, int n_blocks,
                                int block_size, uint64_t master_seed);
```

```text
Same physics as pimc_run, but runs n_replicas fully independent PIMC walkers
(each its own P-bead ring polymer, own equilibration, own block-averaged
sampling) and combines them, parallelized over OpenMP threads when built with
-fopenmp

result.n_blocks = n_blocks * n_replicas (informational total). Falls back
to an all-zero result if n_replicas < 1, usual pimc_run validity conditions
fail for every replica, or on allocation failure
```
