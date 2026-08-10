# Hydrogen Atom

The hydrogen atom is the simplest atomic system with an exact analytical solution. This module provides both a numerical solver for the radial Schrödinger equation and closed‑form wavefunctions and energies.

The code is in `physics/hydrogen.h` and `physics/hydrogen.c`. **All quantities are in SI units** (`hbar`, `mass`, `e_charge`, `eps0` are explicit parameters) - this is one of the few modules that does **not** use atomic units.

## Radial Schrödinger Equation

For a central potential $V(r)$, the radial equation for the reduced wavefunction $u(r) = r R(r)$ is

$$
\left[ -\frac{\hbar^2}{2m} \frac{d^2}{dr^2}
+ \frac{\hbar^2 l(l+1)}{2m r^2}
+ V(r) \right] u(r) = E\, u(r),
$$

With $V(r) = -\dfrac{e^2}{4\pi\varepsilon_0 r}$ for hydrogen.

## Implementation Overview

Three functions are provided:

```c
eigen_t *hydrogen_radial_solve(double *r, int N, int l,
                               double hbar, double mass,
                               double e_charge, double eps0);

double hydrogen_energy_level(int n);

cvector_t *hydrogen_radial_wavefunction(double *r, int N, int n, int l);
```

1. Numerical solver

`hydrogen_radial_solve` builds the finite‑difference Hamiltonian for the given angular momentum ll on the radial grid `r[0..N-1]` (which can be uniform or logarithmic) and returns an `eigen_t` with eigenvalues (energies in Joules) and eigenvectors (the radial functions $u(r)$).

```c
int N = 500;
double r_min = 0.0, r_max = 20.0 * AU_LENGTH;
double *r = linspace(r_min, r_max, N);

eigen_t *sol = hydrogen_radial_solve(r, N, /*l=*/0,
                                     HBAR, M_ELECTRON,
                                     E_CHARGE, EPSILON_0);
// sol->eigenvalues[0] is the 1s energy (~ -2.18e-18 J)
```

The solver uses the generic `central_potential_radial_solve` (from `physics/central_potential.h`) with the Coulomb potential - so it can be easily adapted to other central potentials by changing the potential function.

2. Analytical energy levels

```c
double E_eV = hydrogen_energy_level(1) / E_CHARGE;   // -13.6 eV
double E_J  = hydrogen_energy_level(2);              // -5.44e-19 J
```

3. Analytical radial wavefunction

`hydrogen_radial_wavefunction` returns the **normalised radial** function $R_{nl}(r)$ (not $u(r)$) sampled on the grid $r$. It is built from associated Laguerre polynomials (`core/special/laguerre.c`).

```c
cvector_t *R_10 = hydrogen_radial_wavefunction(r, N, /*n=*/1, /*l=*/0);
// R_10->data[i] = R_{10}(r_i)
```

## Example: Solving and Comparing with Theory

The example `eg_06_hydrogen.c` solves for the lowest s‑states and prints the numerical energies alongside the analytic $E_n=−13.6 \text{eV}/n^2$. It also saves radial probability densities `hydrogen_radial_1.dat`, `hydrogen_radial_2.dat`, etc.

```c
./build/eg_06_hydrogen
```

Typical output (errors < 0.1% for a fine grid):

```
Lowest 5 s-states (l=0):
 n   Numerical E (eV)   Analytical E (eV)   Error (%)
  1   -1.360000e+01      -1.360000e+01       0.00%
  2   -3.400000e+00      -3.400000e+00       0.00%
  3   -1.511111e+00      -1.511111e+00       0.00%
 ...
```

The saved files contain columns: `r (m)`, `R(r)`, `|R|^2`, and `r^2|R|^2` (the radial probability density).

## Validation

The numerical solver is cross‑checked against the analytic energies and wavefunctions in `test_hydrogen.c`. The test confirms:

- Normalisation of the analytic $R_{10}(r) (integral $\int ∣R∣^2 r^2 dr = 1$).
- The 1s energy matches $−13.6$ eV.

## See Also

- [Central Potentials](central_potential.md) - the generic radial solver used underneath.
- [Special Functions](../internals/special.md) - Laguerre polynomials used for analytic wavefunctions.
- [Fine Structure](fine_structure.md/) - adds relativistic corrections to the hydrogen spectrum.
