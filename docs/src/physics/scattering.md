# Scattering Theory

Scattering theory describes the interaction of particles with potentials, including phase shifts, cross sections, and the Born approximation.

## The Schr$\ddot{o}$dinger Equation in Scattering

For a particle incident on a potential $V(r)$, the scattering wave function has asymptotic form:

$$
\psi(\mathbf{r}) \sim e^{ikz} + f(\theta,\phi)\frac{e^{ikr}}{r}
$$

Where $f(\theta,\phi)$ is the scattering amplitude, and the differential cross section is $\frac{d\sigma}{d\Omega} = |f(\theta,\phi)|^2$.

## Partial Wave Expansion

$$
\psi(\mathbf{r}) = \sum_{l=0}^{\infty} i^l(2l+1) R_l(r) P_l(\cos\theta)
$$

The radial functions satisfy:

$$
\left[-\frac{\hbar^2}{2m}\frac{d^2}{dr^2} + \frac{\hbar^2 l(l+1)}{2mr^2} + V(r)\right]R_l(r) = E R_l(r)
$$

With phase shift $\delta_l$ defined by the asymptotic behavior $R_l(r) \sim \frac{1}{kr}\sin(kr - \frac{l\pi}{2} + \delta_l)$, and scattering amplitude $f(\theta) = \frac{1}{k}\sum_l (2l+1)e^{i\delta_l}\sin\delta_l\, P_l(\cos\theta)$.

## Implementation

`physics/scattering.h` implements the phase-shift calculation and the first Born approximation as two independent, unrelated code paths - not a shared radial solver - each taking a [1D-potentials-style](potentials_1d.md) point evaluator directly:

```c
/* Phase shift for partial wave l, via Numerov integration of the radial
   equation + matching to the asymptotic form. Returns delta_l via \atan2,
   in (-\pi, \pi]. */
double phase_shift(int l, double k, potential_fn V, void *params, double r_min,
                   double r_max, int N, double hbar_sq_2m);
```

`k` is the wavenumber ($E = \texttt{hbar\_sq\_2m} \cdot k^2$); `r_min` is a small nonzero starting radius (avoiding the origin, same convention as [Central Potentials](central_potential.md)); `r_max` must extend well beyond the range of `V`; `N` is the number of radial grid points. Internally this uses the [Numerov integrator](../internals/numerov.md), per the header comment, though the exact matching procedure would need `scattering.c` to confirm.

```c
double delta_0 = phase_shift(0, k, V_yukawa, &yukawa_params,
                             1e-4, 50.0, 2000, hbar_sq_2m);
```

### Born Approximation

```c
/* f(\theta) = -(1 / hbar_sq_2m) * (1/q) * \int_0^\inf r * V(r) * \sin(q * r) dr,
    q = 2k * \sin(\theta / 2) */
complex_t born_amplitude(potential_fn V, void *params, double k, double theta,
                         double r_max, int N, double hbar_sq_2m);

// Differential cross section: d\sigma/d\Omega = |f(\theta)|^2
double born_cross_section(complex_t f_theta);
```

This is the standard first-order (weak-potential) approximation to the scattering amplitude - no phase-shift information needed, valid when $V$ is weak enough that the incident wave is only slightly perturbed.

```c
complex_t f = born_amplitude(V_yukawa, &params, k, theta, 50.0, 2000, hbar_sq_2m);
double dsigma_dOmega = born_cross_section(f);
```

> Note: `phase_shift` and `born_amplitude`/`born_cross_section` are independent methods with different validity regimes (partial waves are exact in principle; Born is a weak-potential approximation) - the header gives no combined full partial-wave-sum-to-total-cross-section helper, so converting a set of `phase_shift` results into $\sigma_{\text{tot}} = \frac{4\pi}{k^2}\sum_l(2l+1)\sin^2\delta_l$ would currently be manual.
