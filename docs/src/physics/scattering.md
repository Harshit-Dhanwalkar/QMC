# Scattering Theory

Scattering theory describes the interaction of particles with potentials, including phase shifts, cross sections, and the Born approximation.

## The Schr$\"{o}$dinger Equation in Scattering

For a particle incident on a potential $V(r)$, the scattering wave function has asymptotic form:

$$
\psi(\mathbf{r}) \sim e^{ikz} + f(\theta,\phi)\frac{e^{ikr}}{r}
$$

Where $f(\theta,\phi)$ is the scattering amplitude.

The differential cross section:

$$
\frac{d\sigma}{d\Omega} = |f(\theta,\phi)|^2
$$

## Partial Wave Expansion

The wavefunction is expanded in spherical harmonics:

$$
\psi(\mathbf{r}) = \sum_{l=0}^{\infty} i^l(2l+1) R_l(r) P_l(\cos\theta)
$$

The radial functions satisfy:

$$
\left[-\frac{\hbar^2}{2m}\frac{d^2}{dr^2} + \frac{\hbar^2 l(l+1)}{2mr^2} + V(r)\right]R_l(r) = E R_l(r)
$$

The phase shift $\delta_l$ is defined by the asymptotic behavior:

$$
R_l(r) \sim \frac{1}{kr}\sin(kr - \frac{l\pi}{2} + \delta_l)
$$

The scattering amplitude:

$$
f(\theta) = \frac{1}{k}\sum_{l=0}^{\infty} (2l+1) e^{i\delta_l} \sin\delta_l P_l(\cos\theta)
$$

## Born Approximation

For weak potentials, the scattering amplitude in the first Born approximation:

$$
f(\mathbf{k}', \mathbf{k}) = -\frac{2m}{4\pi\hbar^2}\int e^{-i\mathbf{q}\cdot\mathbf{r}} V(\mathbf{r}) d^3r
$$

Where $\mathbf{q} = \mathbf{k}' - \mathbf{k}$.

### Coulomb Scattering

For the Coulomb potential, the Born approximation gives:

$$
f(\theta) = -\frac{2m e^2}{4\pi\varepsilon_0 \hbar^2}\frac{1}{4k^2\sin^2(\theta/2)}
$$

## Implementation

```c
double born_approximation(const cvector_t *V, const double *r, int N,
                         double theta, double k, double m, double hbar) {
    // q = 2k sin(theta/2)
    double q = 2.0 * k * sin(theta / 2.0);

    // Fourier transform of V(r)
    complex_t integral = c_zero();
    for (int i = 0; i < N; i++) {
        double r_i = r[i];
        complex_t phase = c_exp(c_imag(-q * r_i));
        integral = c_add(integral, c_mul(c_real(V->data[i].re), phase));
    }

    double prefactor = -2.0 * m / (4.0 * M_PI * hbar * hbar);
    return prefactor * integral.re;
}

// Phase shift calculation
double phase_shift(double *r, double *V, int N, int l, double k,
                  double m, double hbar) {
    // Integrate radial Schrödinger equation to large r
    // Match to asymptotic form to extract delta_l
    double R, dRdr;
    integrate_radial(r, V, N, l, k, m, hbar, &R, &dRdr);

    // Extract phase shift from matching at boundary
    double kr = k * r[N-1];
    double delta = atan((kr * R) / (dRdr * r[N-1] - R));
    return delta;
}
```
