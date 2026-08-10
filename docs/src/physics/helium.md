# Helium & Two-Electron Atoms

The Helium atom is a two-electron system ($Z=2$) governed by the non-separable Hamiltonian:

$$
H = -\frac{\hbar^2}{2m} \nabla_1^2 -\frac{\hbar^2}{2m} \nabla_2^2 - \frac{2e^2}{4\pi\varepsilon_0 r_1} - \frac{2e^2}{4\pi\varepsilon_0 r_2} + \frac{e^2}{4\pi\varepsilon_0 \vert{}\mathbf{r}_1 - \mathbf{r}_2\vert{}}
$$

Ground-state energy of a two-electron atom or ion via the classic effective-nuclear-charge variational method (Griffiths Ch. 7; Bransden & Joachain Ch. 5–7). Implemented in `physics/helium.h`.

## The variational ansatz

The trial wavefunction is a product of two hydrogen-like 1s orbitals with a shared _effective_ nuclear charge $Z'$:

$$
\psi(\mathbf{r}_1, \mathbf{r}_2) = \phi(\mathbf{r}_1; Z')\,\phi(\mathbf{r}_2; Z')
$$

> Spatially symmetric, appropriate for the spin-singlet ground state - this models the spatial variational energy only. Electron spin and full antisymmetrization aren't modeled explicitly here (see [Identical Particles](identical.md) for the Slater-determinant machinery that would be needed for that).

All energies are in Hartree atomic units ($\hbar = m_e = e = 4\pi\varepsilon_0 = 1$).

Using known hydrogenic-1s expectation values for the trial charge $Z'$ (valid because $\phi(r;Z')$ is an eigenfunction of the effective one-electron Hamiltonian $H' = -\frac12\nabla^2 - Z'/r$, by the virial theorem):

$$
\langle T\rangle = \frac{Z'^2}{2}\ \text{per electron}, \quad
\left\langle\frac{1}{r}\right\rangle = Z', \quad
\left\langle\frac{1}{r_{12}}\right\rangle = \frac{5}{8}Z'
$$

Gives the closed-form total energy as a function of the trial charge:

$$
E(Z') = Z'^2 - 2ZZ' + \frac{5}{8}Z'
$$

Exactly minimized at $Z'_{\text{opt}} = Z - \frac{5}{16}$, giving
$E_{\text{opt}} = -(Z - \frac{5}{16})^2$. For helium ($Z=2$):
$Z'_{\text{opt}} = 1.6875$, $E_{\text{opt}} = -2.84765625$ Hartree ($\approx -77.5$ eV).

> By the variational theorem, $E_{\text{opt}} \geq E_{\text{true}}$ - the gap to the experimental helium ground-state energy is the known cost of neglecting electron correlation beyond this simple product-orbital ansatz.

## Implementation

```c
// E(Z') in Hartree, for nuclear charge Z (2 for neutral helium)
double helium_variational_energy(double Z_eff, double Z);

// Exact analytic minimizer: Z'_opt = Z - 5/16
double helium_optimal_zeff_analytic(double Z);

// Exact analytic minimum energy: E_opt = -(Z-5/16)^2, in Hartree
double helium_ground_state_energy_analytic(double Z);
```

For validation, or for extending this to a variational form that no longer has a closed-form minimum, there's also a numeric minimizer built on the generic golden-section search (`physics/variational.h`):

```c
/* Numerically minimizes helium_variational_energy over Z' via
   golden_section_minimize, search bounds [max(0.01, 0.1 * Z), 2 * Z].
   If Zeff_opt_out is non-NULL, the optimal Z' found. */
double helium_ground_state_energy_numeric(double Z, double tol,
                                          double *Zeff_opt_out);
```

```c
double Z_eff_found;
double E = helium_ground_state_energy_numeric(2.0, 1e-10, &Z_eff_found);
// E should agree with helium_ground_state_energy_analytic(2.0) to ~tol
```

This numeric/analytic pair is a natural validation check: they should agree to within `tol`, cross-confirming both the closed-form derivation and the golden-section minimizer independently.

## Running the Example

```sh
./build/eg_08_helium
```
