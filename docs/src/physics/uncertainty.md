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

## Numerical Verification

For a given wavefunction, we can compute $\Delta x$ and $\Delta p$ numerically:

```c
double position_uncertainty(const cvector_t *psi, const double *x,
                            int N, double dx) {
    double mean_x = expectation_position(psi, x, N, dx);
    double mean_x2 = expectation_x2(psi, x, N, dx);
    return sqrt(mean_x2 - mean_x*mean_x);
}

double momentum_uncertainty(const cvector_t *psi, const double *x,
                            int N, double dx, double hbar) {
    complex_t mean_p = expectation_momentum(psi, x, N, dx, hbar);
    // Compute <p^2> similarly...
    return sqrt(mean_p2 - c_abs2(mean_p));
}
```

## Minimum Uncertainty States

The Gaussian wavepacket saturates the inequality:

$$
\phi(x)=\frac{1}{(2\pi\sigma^2)^{1/4}} * \exp\left(−\frac{(x−x_0)^2}{4 \sigma^2} + \frac{i p_0 x}{\hbar}\right)
$$

For this state, $\Delta x=\sigma$ and $\Delta p=\hbar/(2\sigma)$, so $\Delta x\Delta p=\hbar/2$.

### Examples

- Infinite square well ground state: $\Delta x\approx 0.18L$, $\Delta p \approx \hbar/L$, $\Delta x \Delta p > \hbar/2$.
- Harmonic oscillator ground state: saturates the bound.
