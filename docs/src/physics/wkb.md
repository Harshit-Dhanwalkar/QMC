# WKB Approximation

The Wentzel-Kramers-Brillouin (WKB) approximation is a semiclassical method for solving the Schrödinger equation when the potential varies slowly compared to the de Broglie wavelength.

## The WKB Wavefunction

The WKB ansatz:

$$
\psi(x) = A \exp\left(\frac{i}{\hbar} S(x)\right)
$$

The semiclassical expansion yields:

$$
\psi(x) \approx \frac{C}{\sqrt{p(x)}} \exp\left(\pm \frac{i}{\hbar}\int^x p(x')dx'\right)
$$

Where the classical momentum is:

$$
p(x) = \sqrt{2m(E - V(x))}
$$

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

The WKB quantization condition:

$$
\int_{x_1}^{x_2} p(x) dx = \left(n + \frac{1}{2}\right)\pi\hbar
$$

Where $x_1$ and $x_2$ are the turning points.

## Implementation

```c
double wkb_integral(const double *x, const double *V, int N,
                    double E, double m, double hbar) {
    double integral = 0.0;
    double p_prev = sqrt(2*m*(E - V[0]));

    for (int i = 1; i < N; i++) {
        double p = sqrt(2*m*(E - V[i]));
        integral += 0.5 * (p + p_prev) * (x[i] - x[i-1]);
        p_prev = p;
    }
    return integral;
}

int wkb_quantization(double *x, double *V, int N, double m, double hbar,
                     int n, double *E_out) {
    // Find energy where integral = (n + 0.5)*pi*hbar
    double E_low = 0.0, E_high = 0.0;
    // Find bounds
    for (int i = 0; i < N; i++) {
        if (E > V[i]) E_high = V[i] + 10.0;
    }
    // Bisection
    for (int iter = 0; iter < 100; iter++) {
        double E_mid = 0.5 * (E_low + E_high);
        double integral = wkb_integral(x, V, N, E_mid, m, hbar);
        double target = (n + 0.5) * M_PI * hbar;
        if (integral < target) E_low = E_mid;
        else E_high = E_mid;
    }
    *E_out = 0.5 * (E_low + E_high);
    return 0;
}
```

## Examples

### Harmonic Osillator

For the harmonic oscillator, $p(x) = \sqrt{2m(E − \frac{1}{2}m\omega^2x^2)}$.

The quantization condition gives:

$$
E_n = \hbar \omega \left(n+\frac{1}{2}\right)
$$

The WKB approximation is exact for the harmonic oscillator.

### Morse Potential

The Morse potential:

$$
V(x) = D_e \left(1 - e^{-a(x-x_0)}\right)^2
$$

has WKB energy levels:

$$
E_n = \hbar \omega \left(n+\frac{1}{2}\right) − \frac{[\hbar \omega(n+1/2)]^2}{4D_e}
$$
