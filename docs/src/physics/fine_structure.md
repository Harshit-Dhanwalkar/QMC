# Fine Structure Corrections

Fine structure corrections account for relativistic effects in single-electron systems, breaking the energy degeneracy with respect to total angular momentum $j$.

## What's actually implemented

The full fine-structure Hamiltonian conceptually has three pieces:

$$H_{\text{fs}} = H_{\text{rel}} + H_{\text{SO}} + H_{\text{Darwin}}$$

`physics/fine_structure.h` implements the **spin-orbit term explicitly** plus a **closed-form total shift** for cross-checking - it does _not_ expose separate functions for $H_{\text{rel}}$ (relativistic kinetic correction) or $H_{\text{Darwin}}$ individually. The background formulas for those two are kept below for context, but treat them as physics reference, not implemented API.

### 1. Relativistic Kinetic Energy Correction (not separately implemented)

$$H_{\text{rel}} = -\frac{p^4}{8m^3c^2} \implies E_{\text{rel}}^{(1)} = -\frac{E_n^2}{2mc^2} \left( \frac{4n}{l + 1/2} - 3 \right)$$

### 2. Spin-Orbit Coupling

$$H_{\text{SO}} = \frac{e^2}{8\pi\varepsilon_0 m^2c^2 r^3} \mathbf{L} \cdot \mathbf{S}$$

### 3. Darwin Term ($l=0$)

$$H_{\text{Darwin}} = \frac{\pi\hbar^2 e^2}{2m^2c^2 (4\pi\varepsilon_0)} \delta^3(\mathbf{r})$$

## Implementation

Quantum numbers here follow `angular.c`'s **doubled** convention throughout - `j_2 = 2j`, `M_2 = 2M` - so half-integer $j$ (from coupling orbital $l$ with spin-$\frac12$) stays an integer in code. For an electron, $j = l \pm \frac12 \Rightarrow j_2 = 2l \pm 1$.

```c
// <L.S>/\hbar^2 = (1/2) * (j(j+1) - l(l+1) - s(s+1)), s=1/2 fixed. Closed form
double spin_orbit_ls_expect(int l, int j_2);
```

$$
\frac{\langle \mathbf{L}\cdot\mathbf{S}\rangle}{\hbar^2} = \frac{1}{2}\left[j(j+1) - l(l+1) - s(s+1)\right], \quad s = \tfrac12
$$

```c
// <1/r^3>_{nl} for hydrogen = 1/(a0^3 n^3 l(l+1/2)(l+1)). Returns 0 for l=0
double hydrogen_expect_inv_r3(int n, int l, double hbar, double mass,
                              double e_charge, double eps0);
```

```c
/* Perturbative spin-orbit energy shift. Returns 0 for l=0.
   spin-orbit contribution only - total physical fine structure
   also needs relativistic kinetic correction and (for l=0) Darwin
   term, neither of which this function includes. */
double hydrogen_spin_orbit_energy(int n, int l, int j_2, double hbar,
                                  double mass, double e_charge, double eps0,
                                  double c);
```

$$
\delta E_{\text{SO}} = \frac{e^2}{8\pi\varepsilon_0 m^2c^2}\left\langle\frac{1}{r^3}\right\rangle_{nl} \hbar^2\,\frac{\langle\mathbf{L}\cdot\mathbf{S}\rangle}{\hbar^2}
$$

### Combined Energy Shift (closed form, hydrogen only)

```c
// dE_fs = E_n * (\alpha^2 / n^2) * (n / (j+1/2) - 3/4)
double hydrogen_fine_structure_shift(int n, int j_2, double hbar, double mass,
                                     double e_charge, double eps0, double c);
```

$$
E_{\text{fs}} = \frac{E_n\, \alpha^2}{n^2}\left(\frac{n}{j+1/2} - \frac{3}{4}\right)
$$

Where $\alpha \approx 1/137.036$. Note this signature has **no `Z` parameter** - it's the pure-hydrogen ($Z=1$) closed form, algebraically equivalent to $E_n\alpha^2 Z^2/n \cdot (1/(j+\frac12) - 3/(4n))$ with $Z=1$ substituted in. If you need the $Z$-generalized hydrogenic version, it isn't exposed here.

### Independent Validation

```c
/* Recomputes <L.S> from scratch via couple_states()'s Clebsch-Gordan coefficients and l_plus_op/l_minus_op, as a cross-check against closed-form spin_orbit_ls_expect() above. */
double spin_orbit_ls_expect_from_coupling(int l, int j_2, int M_2);
```

Evaluates $\mathbf{L}\cdot\mathbf{S} = L_zS_z + \frac12(L_+S_- + L_-S_+)$ term-by-term on the uncoupled product basis using [angular momentum coupling](angular_momentum.md) machinery, giving an independent numerical check that the ladder-operator and Clebsch-Gordan implementations agree with the closed-form $j(j+1)-l(l+1)-s(s+1)$ result. Returns `NAN` for invalid `(l, j_2, M_2)`.
