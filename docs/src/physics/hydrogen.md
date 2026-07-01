# Hydrogen Atom

The hydrogen atom is the quantum mechanical solution of a Coulomb potential $V(r) = -e^2/(4\pi\varepsilon_0 r)$. It's the only atomic system with exact analytical solutions.

## The Schr$\"{o}$dinger Equation

The time-independent Schr$\"{o}$dinger equation in spherical coordinates:

$$
\left[-\frac{\hbar^2}{2m}\nabla^2 - \frac{e^2}{4\pi\varepsilon_0 r}\right]\psi = E\psi
$$

The wavefunction separates:

$$
\psi_{nlm}(r, \theta, \phi) = R_{nl}(r) Y_{lm}(\theta, \phi)
$$

where $Y_{lm}$ are spherical harmonics and $R_{nl}$ are radial functions.

### Radial Equation

$$
-\frac{\hbar^2}{2m}\frac{d^2R}{dr^2} - \frac{\hbar^2}{mr}\frac{dR}{dr} + \left[\frac{\hbar^2 l(l+1)}{2mr^2} - \frac{e^2}{4\pi\varepsilon_0 r}\right]R = ER
$$

The bound state energies are:

$$
E_n = -\frac{m e^4}{2(4\pi\varepsilon_0)^2\hbar^2}\frac{1}{n^2} = -13.6\,\text{eV}\frac{1}{n^2}
$$

### Radial Wavefunctions

For $l = 0$ (s-states):

$$
R*{n0}(r) = 2\left(\frac{1}{na_0}\right)^{3/2} \sqrt{\frac{1}{n^2}}
L*{n-1}^1(2r/na_0) e^{-r/na_0}
$$

For $l = 1$ (p-states):

$$
R\_{21}(r) = \frac{1}{2\sqrt{6}}\left(\frac{1}{a_0}\right)^{3/2}
\frac{r}{a_0} e^{-r/2a_0}
$$

## Numerical Solution

The radial equation is solved in `src/physics/hydrogen.c`:

```c
eigen_t *hydrogen_radial_solve(double *r, int N, int l,
                              double hbar, double m, double e, double eps0) {
    // Build Hamiltonian matrix for radial equation
    cmatrix_t *H = cmatrix_alloc(N, N);

    double coeff = hbar*hbar / (2*m*dr*dr);

    for (int i = 0; i < N; i++) {
        double V_coulomb = -e*e / (4*M_PI*eps0*r[i]);
        double V_centrifugal = hbar*hbar * l*(l+1) / (2*m*r[i]*r[i]);

        CMAT(H, i, i) = c_real(2*coeff + V_coulomb + V_centrifugal);
        if (i > 0) CMAT(H, i, i-1) = c_real(-coeff);
        if (i < N-1) CMAT(H, i, i+1) = c_real(-coeff);
    }

    return cmatrix_eigh(H);
}
```

### Running the Example

```sh
./build/eg_03_hydrogen
```

### Output

// TODO

The radial probability density $r^2∣R_{nl}(r)∣^2$ is saved to `hydrogen_radial_*.dat`.
