# Path Integral Monte Carlo (PIMC)

Path Integral Monte Carlo samples the thermal density matrix of a quantum system via the Feynman path integral. At low temperatures (large imaginary time $\beta$), the thermodynamic energy approaches the ground‑state energy without relying on a trial wavefunction. For two‑electron atoms/ions, the method treats the electrons as distinguishable particles (Boltzmannons) and uses the Kelbg‑regularised Coulomb potential to avoid the path‑collapse catastrophe.

Implemented in `physics/pimc.h` and `physics/pimc.c`. All energies in Hartree atomic units.

## Path Integral Representation

For two electrons in a central potential, the thermal density matrix is written as a ring polymer of $P$ beads (imaginary‑time slices). The action for each bead is

$$
S = \tau \sum_{k=1}^P \left[ \frac{(\mathbf{r}_{1,k} - \mathbf{r}_{1,k+1})^2}{2\tau^2}
+ \frac{(\mathbf{r}_{2,k} - \mathbf{r}_{2,k+1})^2}{2\tau^2}
+ V_{\rm Kelbg}(\mathbf{r}_{1,k}, \mathbf{r}_{2,k}) \right],
$$

Where $\tau = \beta/P$ and $V_{\rm Kelbg}$ is the Kelbg pair potential.

### Kelbg Potential

The bare Coulomb potential has a $1/r$ singularity that causes numerical instability in path‑integral sampling. The Kelbg regularisation replaces it with

$$
V_{\rm Kelbg}(r, \lambda, q) =
\frac{q}{r}\left(1 - e^{-r^2/\lambda^2}\right)
+ \frac{q\sqrt{\pi}}{\lambda}\, \operatorname{erfc}\!\left(\frac{r}{\lambda}\right),
$$

Where $q = q_a q_b$ and $\lambda = \sqrt{\tau/\mu}$ (with $\mu$ the reduced mass). This potential is finite at $r=0$ and reduces to the bare Coulomb as $\tau \to 0$.

The function `kelbg_potential()` implements this; `kelbg_energy_correction()` provides the $\tau$-dependent part needed for the thermodynamic energy estimator.

### Bisection Sampling

To sample the path, we use **bisection moves** (Ceperley 1995). A segment of length $2^L$ beads (from bead $i$ to $i+2^L$) is resampled from the exact free‑particle bridge distribution, and the whole segment is accepted or rejected based on the potential‑action difference. This is implemented in `pimc_bisection_move()`.

### Energy Estimator

The thermodynamic energy estimator is

$$
E = \frac{3N}{2\tau} - \frac{1}{P}\sum_{k=1}^P
\left[
\frac{(\Delta \mathbf{r}_{1,k})^2 + (\Delta \mathbf{r}_{2,k})^2}{2\tau^2}
+ V_{\rm Kelbg}^{\rm corr}(\mathbf{r}_{1,k}, \mathbf{r}_{2,k})
\right],
$$

With $N=2$ electrons and $V_{\rm Kelbg}^{\rm corr}$ given by `kelbg_energy_correction()`. The estimator is evaluated by `pimc_energy_estimator()`.

## API Overview

```c
typedef struct {
    int P;
    double (*r1)[3];
    double (*r2)[3];
} pimc_walker_t;

typedef struct {
    double energy;
    double error;
    double acceptance_rate;
    int n_blocks;
} pimc_result_t;

pimc_result_t pimc_run(double Z, int P, double tau, int level,
                       int n_equilibration, int n_blocks, int block_size,
                       uint64_t seed);
```

## Example: Helium Ground State

```c
#include "physics/pimc.h"

double Z = 2.0;
double beta = 8.0;       // low temperature
int P = 512;             // must be power of 2
double tau = beta / P;
int level = 4;           // bisection segment length = 16 beads

pimc_result_t r = pimc_run(Z, P, tau, level,
                           400, 30, 250,   // equilibration, blocks, block size
                           20260802ULL);
printf("PIMC energy: %f +- %f Hartree\n", r.energy, r.error);
// Typical output: ~ -2.85 ± 0.05 Hartree (depending on P and statistics)
```

## Validation

Tests (`test_pimc.c`) check:

- The Kelbg potential and energy correction against closed‑form values.
- The potential is finite at r=0r=0 and reduces to the Coulomb limit.
- Invalid inputs (non‑power‑of‑two PP, out‑of‑range level) are rejected.
- A full PIMC run at P=512P=512, τ=0.015625τ=0.015625 gives an energy within ~0.4 Hartree of the exact ground state (accepting that the thermodynamic estimator variance grows with PP; larger statistics or a virial estimator would improve this).

## See Also

- [Variational Monte Carlo](vmc.md) - uses a trial wavefunction.
- [Diffusion Monte Carlo](dmc.md) - importance‑sampled projection to the ground state.
