# Diffusion Monte Carlo (DMC)

Diffusion Monte Carlo projects the ground state from a trial wavefunction by evolving in imaginary time. With importance sampling (using the trial wavefunction as a guide), the method can yield the exact ground‑state energy (within statistical and finite‑time‑step errors) without the fixed‑node approximation for systems with no nodes - such as the helium ground state.

Implemented in `physics/dmc.h` and `physics/dmc.c`. All energies in Hartree atomic units.

## Importance‑Sampled DMC

The imaginary‑time Schrödinger equation with importance sampling is transformed into a drift‑diffusion process for walkers:

$$
\mathbf{r}' = \mathbf{r} + \tau \, \mathbf{v}(\mathbf{r}) + \chi, \qquad \chi \sim \mathcal{N}(0, \tau\, \mathbf{I}_3),
$$

Where the drift velocity is

$$
\mathbf{v}(\mathbf{r}) = \nabla \ln |\Psi_T(\mathbf{r})|^2 = 2\,\nabla \ln \Psi_T(\mathbf{r}).
$$

For the Slater‑Jastrow trial wavefunction (same as in VMC), the drift velocity is computed analytically in `dmc_drift_velocity()`.

### Move Acceptance

Each single‑electron move is accepted with the Metropolis‑corrected probability that accounts for the drift‑diffusion Green's function:

$$
A = \min\!\left(1,\;
\frac{|\Psi_T(\mathbf{r}')|^2\, G(\mathbf{r}\leftarrow \mathbf{r}')}
     {|\Psi_T(\mathbf{r})|^2\, G(\mathbf{r}'\leftarrow \mathbf{r})}
\right),
$$

Where $G$ is the drift‑diffusion propagator.

### Branching

After moving both electrons, each walker is assigned a multiplicity

$$
m = \text{round}\!\left[ \exp\!\left( -\tau \left( \frac{E_L(\text{old}) + E_L(\text{new})}{2} - E_T \right) \right) \right],
$$

Where $E_T$ is a reference energy adjusted to keep the population stable. Walkers with $m=0$ are killed; those with $m>1$ are replicated.

Population is controlled by a feedback loop on $E_T$ and by uniform subsampling when the population exceeds a maximum cap.

### Estimators

Two estimators for the ground‑state energy are available:

- **Mixed estimator**: the average local energy over the walker population.
- **Growth estimator**: the time‑averaged $E_T$ history (also block‑averaged).

The mixed estimator has smaller variance and is the standard DMC result.

## API Overview

```c
typedef struct {
    vmc_walker_t *data;
    int count;
    int capacity;
} dmc_population_t;

typedef struct {
    double energy_mixed;     // block‑averaged mixed estimator
    double error_mixed;
    double energy_growth;    // block‑averaged growth estimator
    double error_growth;
    int n_blocks;
    double mean_population;
    double acceptance_rate;
} dmc_result_t;

dmc_result_t dmc_run(double Z, double Zeff, double b,
                     int target_population, int max_population,
                     double tau, int n_equilibration,
                     int n_blocks, int block_size,
                     uint64_t seed);
```

## Example: Helium Ground State

```c
#include "physics/dmc.h"

double Z = 2.0, Zeff = 2.0, b = 0.15;
dmc_result_t r = dmc_run(Z, Zeff, b,
                         500, 1500,    // target and max population
                         0.01,         // time step
                         1000,         // equilibration generations
                         30, 200,      // 30 blocks × 200 generations
                         5678ULL);
printf("Mixed estimator: %f +- %f Hartree\n", r.energy_mixed, r.error_mixed);
// Typical output: ~ -2.891 ± 0.005 Hartree
```

Full run is shown in `eg_29_dmc_helium.c`.

## Validation

Tests (`test_dmc.c`) verify:

- Deterministic drift‑velocity values against analytical derivatives.
- The variational theorem: DMC energy is below the VMC energy for the same trial wavefunction.
- The mixed estimator is close to the exact non‑relativistic helium energy ($−2.9037$ Hartree) within a tolerance that accounts for finite time step and population bias.
- Population control keeps the average walker count near the target.

## See Also

- [Variational Monte Carlo](vmc.md) - the sampling method used to initialise and guide DMC.
- [Path Integral Monte Carlo](pimc.md) - a different stochastic approach without importance sampling.
