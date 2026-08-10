# Variational Monte Carlo (VMC)

Variational Monte Carlo uses stochastic sampling to evaluate the expectation value of the Hamiltonian for a parameterised trial wavefunction. For two‑electron atoms/ions, the method provides a variational upper bound to the ground‑state energy and can optimise correlation parameters (e.g., the Jastrow parameter $b$) to recover a significant fraction of the correlation energy.

Implemented in `physics/vmc.h` and `physics/vmc.c`. All energies are in Hartree atomic units ($\hbar = m_e = e = 4 * \pi * \varepsilon_0 = 1$).

## Trial Wavefunction

The Slater‑Jastrow ansatz for a two‑electron atom/ion with nuclear charge $Z$ is

$$
\Psi_T(\mathbf{r}_1, \mathbf{r}_2) = e^{-Z_{\rm eff}(r_1 + r_2)} \,
\exp\!\left( \frac{r_{12}}{2(1 + b\,r_{12})} \right),
$$

where:

- $r_i = |\mathbf{r}_i|$
- $r_{12} = |\mathbf{r}_1 - \mathbf{r}_2|$
- $Z_{\rm eff}$ is a variational orbital exponent (usually set to the bare nuclear charge or optimised separately)
- $b$ is the Jastrow parameter controlling the strength of the electron‑electron cusp.

The prefactor $1/2$ in the Jastrow exponent is fixed by the electron‑electron cusp condition; only $b$ is varied.

## Local Energy

For a given configuration, the local energy is defined as

$$
E_L(\mathbf{R}) = \frac{\hat{H}\Psi_T(\mathbf{R})}{\Psi_T(\mathbf{R})}.
$$

For the above wavefunction, the local energy is computed analytically (not by finite differences) and is implemented in `vmc_local_energy()`. The expression is

$$
\begin{aligned}
E_L ={} & -Z_{\rm eff}^2
+ (Z_{\rm eff} - Z)\left(\frac{1}{r_1}+\frac{1}{r_2}\right) \\
& + \frac{b}{(1+b s)^3} - \frac{1}{s(1+b s)^2}
- \frac{1}{4(1+b s)^4} \\
& + \frac{Z_{\rm eff}}{2(1+b s)^2}
\left( \hat{\mathbf{r}}_1 - \hat{\mathbf{r}}_2 \right)\cdot \hat{\mathbf{r}}_{12}
+ \frac{1}{s},
\end{aligned}
$$

With $s = r_{12}$.

## Monte Carlo Sampling

The configuration space is sampled according to $|\Psi_T|^2$ using the Metropolis algorithm.

- A **walker** stores the positions of both electrons: `vmc_walker_t`.
- Each **sweep** attempts to move one electron at a time with a uniform random displacement in a cube of side $2 \times \text{step\_size}$.
- The move is accepted with probability $\min(1, |\Psi_T(\text{new})|^2 / |\Psi_T(\text{old})|^2)$.

After equilibration, the local energy is accumulated and block‑averaged to estimate $\langle E_L \rangle$ and its standard error.

## API Overview

```c
typedef struct {
    double r1[3];   // electron 1 position
    double r2[3];   // electron 2 position
} vmc_walker_t;

typedef struct {
    double mean;             // <E_L>
    double error;            // standard error
    double variance;         // sample variance
    int n_samples;
    double acceptance_rate1;
    double acceptance_rate2;
} vmc_result_t;

double vmc_trial_wavefunction(const vmc_walker_t *w, double Zeff, double b);
double vmc_local_energy(const vmc_walker_t *w, double Z, double Zeff, double b);
void vmc_walker_init(vmc_walker_t *w, rng_state_t *rng, double Zeff);
int vmc_metropolis_move_electron(vmc_walker_t *w, int which, double Zeff,
                                 double b, double step_size, rng_state_t *rng);
void vmc_metropolis_sweep(vmc_walker_t *w, double Zeff, double b,
                          double step_size1, double step_size2,
                          rng_state_t *rng, int *accepted1, int *accepted2);
vmc_result_t vmc_run(double Z, double Zeff, double b, int n_equilibration,
                     int n_samples, int block_size, double step_size1,
                     double step_size2, uint64_t seed);
double vmc_optimize_b(double Z, double Zeff, double b_min, double b_max,
                      int n_equilibration, int n_samples, double step_size1,
                      double step_size2, uint64_t seed, double tol,
                      double *b_opt_out);
```

## Example: Helium Ground State

```c
#include "physics/vmc.h"
#include "physics/helium.h"

double Z = 2.0, Zeff = 2.0;
double b_opt;
double E_opt = vmc_optimize_b(Z, Zeff, 0.0, 0.6, 1000, 50000,
                              0.9, 0.9, 4242ULL, 1e-3, &b_opt);
printf("b_opt = %f, E = %f Hartree\n", b_opt, E_opt);
// Typical output: b_opt ≈ 0.18, E ≈ -2.878 Hartree
```

The complete sweep over _b_ is shown in `eg_28_vmc_helium.c`.

## Validation

The deterministic local‑energy values are checked against independent analytical references (cross‑verified with finite‑difference Laplacian) in `test_vmc.c`. The VMC energy satisfies the variational theorem:

$$
E_{VMC} \ge E{exact} (within statistical error),
$$

and the Jastrow factor lowers the energy below the simple product‑orbital result ($−2.84765625$ Hartree) toward the exact value ($−2.9037$ Hartree).

## See Also

- [Helium & Two‑Electron Atoms](helium.md) - the simpler variational approach without Monte Carlo.
- [Diffusion Monte Carlo](dmc.md) - projects beyond the variational limit.
- [Path Integral Monte Carlo](pimc.md) - a different, trial‑wavefunction‑free approach.
