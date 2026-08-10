# Uncertainty Principle

The Heisenberg uncertainty principle states that certain pairs of observables cannot be simultaneously known with arbitrary precision.

## Position-Momentum Uncertainty

$$
\Delta x \cdot \Delta p \geq \frac{\hbar}{2}
$$

Where $\Delta x = \sqrt{\langle x^2 \rangle - \langle x \rangle^2}$ and similarly for $p$.

## Energy-Time Uncertainty

$$
\Delta E \cdot \Delta t \geq \frac{\hbar}{2}
$$

This governs the lifetime of quantum states and the natural linewidth.

## Implementation

`physics/uncertainty.h` computes everything in one call against a
[`wavefunction_t`](wavefunctions.md), rather than separate position/momentum
functions:

```c
typedef struct {
  double mean_x;
  double mean_p;
  double delta_x;
  double delta_p;
  double product; // \Delta x * \Delta p
} uncertainty_t;

uncertainty_t compute_uncertainties(const wavefunction_t *wf);
```

```c
wavefunction_t *wf = /* build and normalize */;
uncertainty_t u = compute_uncertainties(wf);

printf("<x>=%g  <p>=%g  dx=%g  dp=%g  dx*dp=%g\n",
       u.mean_x, u.mean_p, u.delta_x, u.delta_p, u.product);
```

Internally this is presumably built from `wavefunction_expect_x`,
`wavefunction_expect_x2`, `wavefunction_expect_p`, and `wavefunction_expect_p2`
(see [Wave Functions](wavefunctions.md)), though the exact implementation
would need `uncertainty.c` to confirm.

### Energy expectation

The same header also provides the energy expectation value for a given
potential and mass, independent of the uncertainty struct above:

```c
double compute_energy_expectation(const wavefunction_t *wf, potential_fn V,
                                  void *params, double mass);
```

Takes a `potential_fn` in the same point-evaluator form as
[1D Potentials](potentials_1d.md) (e.g. `V_harmonic`), plus its `params`, and
the particle mass explicitly - since `wavefunction_t` itself carries no mass
information.

## Minimum Uncertainty States

The Gaussian wavepacket saturates the inequality:

$$
\phi(x)=\frac{1}{(2\pi\sigma^2)^{1/4}} * \exp\left(−\frac{(x−x_0)^2}{4 \sigma^2} + \frac{i p_0 x}{\hbar}\right)
$$

For this state, $\Delta x=\sigma$ and $\Delta p=\hbar/(2\sigma)$, so $\Delta x\Delta p=\hbar/2$.

### Examples

- Infinite square well ground state: $\Delta x\approx 0.18L$, $\Delta p \approx \hbar/L$, $\Delta x \Delta p > \hbar/2$.
- Harmonic oscillator ground state: saturates the bound.
