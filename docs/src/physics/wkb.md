# WKB Approximation

The Wentzel-Kramers-Brillouin (WKB) approximation is a semiclassical method for solving the Schr$\ddot{o}$dinger equation when the potential varies slowly compared to the de Broglie wavelength.

## The WKB Wavefunction

The WKB ansatz:

$$
\psi(x) = A \exp\left(\frac{i}{\hbar} S(x)\right)
$$

The semiclassical expansion yields:

$$
\psi(x) \approx \frac{C}{\sqrt{p(x)}} \exp\left(\pm \frac{i}{\hbar}\int^x p(x')dx'\right)
$$

Where the classical momentum is $p(x) = \sqrt{2m(E - V(x))}$.

## Regions of Validity

The WKB approximation is valid when:

$$
\left|\frac{d}{dx}\left(\frac{\hbar}{p(x)}\right)\right| \ll 1
$$

This fails near classical turning points where $E \approx V(x)$.

## Connection Formulas

Near a turning point, the WKB solutions must be matched using Airy functions. For a right turning point at $x = a$:

**For $x < a$** (classically allowed):

$$
\psi(x) = \frac{C}{\sqrt{p(x)}} \cos\left(\frac{1}{\hbar}\int_x^a p(x')dx' - \frac{\pi}{4}\right)
$$

**For $x > a$** (classically forbidden):

$$
\psi(x) = \frac{C}{2\sqrt{|p(x)|}} \exp\left(-\frac{1}{\hbar}\int_a^x |p(x')|dx'\right)
$$

## Bohr-Sommerfeld Quantization

$$
\int_{x_1}^{x_2} p(x) dx = \left(n + \frac{1}{2}\right)\pi\hbar
$$

## Implementation

`physics/wkb.h` takes a [1D-potentials-style](potentials_1d.md) point evaluator (`potential_fn V(x, params)`) directly, rather than a pre-filled array - the quantization/transmission integrals are computed internally against the potential function, not a discretized `V[]`:

```c
double wkb_quantization(potential_fn V, void *params, double m, double hbar,
                        int n, double E_min, double E_max, double x_min,
                        double x_max, int N, double tol);
```

Solves $\int_{x_1}^{x_2}\sqrt{2m(E-V(x))}\,dx = (n+\frac12)\pi\hbar$ for the energy of level `n`, bracketed on `[E_min, E_max]`, integrated over `[x_min, x_max]` with `N` points, to tolerance `tol`. Turning points $x_1, x_2$ are presumably located internally from where $E = V(x)$, rather than being caller-supplied.

```c
double wkb_transmission(potential_fn V, void *params, double m, double hbar,
                        double E, double x1, double x2, int N);
```

$$
T = \exp\left(-\frac{2}{\hbar}\int_{x_1}^{x_2}\sqrt{2m(V(x)-E)}\,dx\right)
$$

Here the turning points `x1, x2` **are** caller-supplied - this is the tunneling-through-a-barrier case, not a bound-state case, so there's no ambiguity about which crossing points to use.

### Closed-form wrappers

Three convenience functions wrap the general solvers above for cases with
known analytic answers - useful both for quick lookups and as validation
targets for the general integrators:

```c
// E_n = \hbar * \omega * (n+1/2), for V(x) = (1/2) * \omega^2 * x^2 (m=1)
double wkb_energy_harmonic(int n, double hbar, double omega);

/* Full round-trip action integral 2 * \int_{-x_turn}^{x_turn} \sqrt(2(E-V(x))) dx
   for the harmonic oscillator (m=1). Quantized levels satisfy
   action / (2 * \pi * \hbar) = n + 1/2. */
double wkb_action_integral_harmonic(double E, double hbar, double omega);

/* Closed-form rectangular-barrier tunneling: T = \exp(-2 * kappa * width),
   \kappa = \sqrt((V0 - E) / hbar_sq_2m).
   Returns 1.0 if V0 <= E.
*/
double wkb_tunneling_rectangular(double E, double V0, double width,
                                 double hbar_sq_2m);
```

`wkb_energy_harmonic` reproduces the famous fact that **WKB is exact for the harmonic oscillator** - the quantization condition gives $E_n = \hbar\omega(n+\frac12)$ with no approximation error, unlike a generic potential where WKB is only asymptotically accurate.

## Examples

### Harmonic Oscillator

For the harmonic oscillator, $p(x) = \sqrt{2m(E - \frac12 m\omega^2 x^2)}$,
and the quantization condition gives the exact spectrum above.

### Morse Potential

The Morse potential $V(x) = D_e(1-e^{-a(x-x_0)})^2$ has WKB energy levels:

$$
E_n = \hbar\omega\left(n+\frac12\right) - \frac{[\hbar\omega(n+1/2)]^2}{4D_e}
$$

> No dedicated Morse wrapper exists in `wkb.h` - use the general `wkb_quantization` with `V_morse` (see [1D Potentials](potentials_1d.md)) to reproduce this formula numerically.
