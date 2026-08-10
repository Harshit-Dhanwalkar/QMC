# Variational Method

The Rayleigh-Ritz variational method provides an upper bound on the true ground state energy $E_0$.

The generic infrastructure behind [Helium](helium.md)'s numeric cross-check and any other single-parameter trial-wavefunction minimization. Implemented in `physics/variational.h`.

## Rayleigh Quotient

For any normalized trial wave function $\psi_{\text{trial}}(\mathbf{r}; \boldsymbol{\alpha})$ containing adjustable parameters $\boldsymbol{\alpha}$:

$$
E(\boldsymbol{\alpha}) = \frac{\langle \psi_{\text{trial}}(\boldsymbol{\alpha}) \vert{} H \vert{} \psi_{\text{trial}}(\boldsymbol{\alpha}) \rangle}{\langle \psi_{\text{trial}}(\boldsymbol{\alpha}) \vert{} \psi_{\text{trial}}(\boldsymbol{\alpha}) \rangle} \ge E_0
$$

The optimal parameters are found by minimizing $E(\boldsymbol{\alpha})$:

$$
\frac{\partial E(\boldsymbol{\alpha})}{\partial \alpha_k} = 0
$$

## Linear Variational Method

Expanding trial function in $M$ basis functions:

$$\psi_{\text{trial}} = \sum_{j=1}^M c_j \phi_j$$

reduces parameter search to the generalized matrix eigenvalue problem:

$$\mathbf{H} \mathbf{c} = E \mathbf{S} \mathbf{c}$$

Where

- $H_{ij} = \langle \phi_i \vert{} H \vert{} \phi_j \rangle$
- $S_{ij} = \langle \phi_i \vert{} \phi_j \rangle$

## Generic 1D minimizer

```c
// Minimizes an arbitrary scalar function f(x, params) over [a,b].
double golden_section_minimize(double a, double b,
                               double (*f)(double x, void *params),
                               void *params, double tol);
```

Golden-section search - derivative-free, guaranteed convergence for a unimodal function on a bounded interval, at the cost of being slower than a gradient-based method. This is what [`helium_ground_state_energy_numeric`](helium.md) calls under the hood to minimize over the effective charge $Z'$, and is generic enough to reuse for any other single-parameter variational problem.

## Trial-wavefunction energy

```c
// Expectation value of H for a trial wavefunction against potential V.
double variational_energy(const wavefunction_t *wf, potential_fn V,
                          void *params, double mass);
```

Effectively the same computation as [`compute_energy_expectation`](uncertainty.md) in `uncertainty.h` - both take a [`wavefunction_t`](wavefunctions.md), a [point-evaluator potential](potentials_1d.md), and a mass, and return $\langle\psi|H|\psi\rangle$. Whether one wraps the other or they're independent implementations of the same formula would need `variational.c` and `uncertainty.c` to confirm; worth resolving so there's a single source of truth for this calculation rather than two.

## Minimizing over a trial-function family

```c
typedef struct {
  void (*trial_func)(double alpha, wavefunction_t *wf);
  wavefunction_t *wf;
  potential_fn V;
  void *params;
  double mass;
} variational_closure_t;

double variational_minimize(double alpha_min, double alpha_max,
                            void (*trial_func)(double alpha,
                                               wavefunction_t *wf),
                            wavefunction_t *wf, potential_fn V, void *params,
                            double mass, double tol);
```

> `trial_func` fills `wf` in place given a trial parameter `alpha` (e.g. a Gaussian width, an effective charge, a variational exponent); `variational_minimize` sweeps `alpha` over `[alpha_min, alpha_max]` via `golden_section_minimize`, evaluating `variational_energy` at each trial point, and returns the minimized energy. `variational_closure_t` looks like the bundle of state this needs to pass through to `golden_section_minimize`'s single-argument `f(x, params)` signature - worth confirming it's actually used that way in `variational.c` rather than being vestigial.

```c
void gaussian_trial(double alpha, wavefunction_t *wf) {
    for (int i = 0; i < wf->n; i++)
        wf->psi->data[i] = c_real(exp(-alpha * wf->x[i] * wf->x[i]));
    wavefunction_normalize(wf);
}

double E_min = variational_minimize(0.1, 5.0, gaussian_trial, wf,
                                    V_harmonic, &omega, /*mass=*/1.0, 1e-10);
```

This general parametrized-trial-function path is the natural next step beyond [Helium](helium.md)'s single closed-form-checkable case - useful anywhere a trial wavefunction's shape (not just one scalar charge) needs tuning against a potential with no exact solution.
