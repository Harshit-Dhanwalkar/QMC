# Angular Momentum and Spin

Quantum mechanical treatment of angular momentum including orbital angular momentum, spin, and their coupling.
This module provides tools for orbital angular momentum, spin-1/2, Clebsch-Gordan coupling, and Gaunt coefficients.

Implemented in `physics/angular.h` and `physics/angular.c`. All quantities are dimensionless (angular momenta in units of $\hbar$).

## Orbital Angular Momentum

The operators $L^2$, $L_z$, $L_x$, $L_y$ are represented as $(2l+1)\times(2l+1)$ matrices in the basis $\{|l,m\rangle\}$, with $m = -l, -l+1, \dots, l$.

The orbital angular momentum operator:

$$
\mathbf{L} = \mathbf{r} \times \mathbf{p}
$$

In spherical coordinates:

$$
L_z = -i\hbar\frac{\partial}{\partial\phi}
$$

$$
L^2 = -\hbar^2\left[\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right) + \frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\phi^2}\right]
$$

### Eigenvalues

$$
L^2|l,m\rangle = \hbar^2 l(l+1)|l,m\rangle, \qquad L_z|l,m\rangle = \hbar m|l,m\rangle
$$

Where $l = 0, 1, 2, \ldots$ and $m = -l, -l+1, \ldots, l$.

### Matrix representations

`physics/angular.h` builds the $(2l+1)\times(2l+1)$ matrix representation of each operator for a fixed $l$, in the $|l,m\rangle$ basis:

```c
cmatrix_t *l2_matrix(int l);
cmatrix_t *lz_matrix(int l);
cmatrix_t *lx_matrix(int l);
cmatrix_t *ly_matrix(int l);
```

Built from the ladder operators:

```c
// <l,m+1|L_+|l,m>
complex_t l_plus_op(int l, int m, int m_prime);
// <l,m-1|L_-|l,m>
complex_t l_minus_op(int l, int m, int m_prime);
```

Where $L_\pm = L_x \pm iL_y$, so $L_x = \frac12(L_++L_-)$ and $L_y = \frac{1}{2i}(L_+-L_-)$ - which is presumably how `lx_matrix` and `ly_matrix` are actually built internally.

Example - building $L_z$ for $l=1$:

```c
cmatrix_t *Lz = lz_matrix(1);   // 3x3 matrix with diag(-1,0,1)
```

## Spherical Harmonics

Spherical harmonics are implemented in core/special/special.h as

```c
complex_t spherical_harmonic(int l, int m, double theta, double phi);
```

They are used in the Gaunt coefficient routine (see below).

### Spin-1/2

The Pauli matrices are stored as flat arrays (row‑major) of `complex_t`:

```c
extern complex_t sigma_x[4];
extern complex_t sigma_y[4];
extern complex_t sigma_z[4];
```

Spinors are plain `cvector_t` of length 2. In‑place operations apply a Pauli matrix:

```c
void spin_sigma_x(cvector_t *spinor);
void spin_sigma_y(cvector_t *spinor);
void spin_sigma_z(cvector_t *spinor);
```

```c
cvector_t *spin = cvector_alloc(2);
spin->data[0] = c_real(1.0);   // |up>
spin->data[1] = c_real(0.0);

spin_sigma_z(spin);             // unchanged (eigenvalue +1)
spin_sigma_x(spin);             // now |down> (swap components)
```

To compute expectation values, use `cvector_dot` or `cvector_expect` from the [core linear algebra module](../internals/linalg.md).

### Coupling of Angular Momenta

Two angular momenta $j_1$ and $j_2$ are coupled to total $J=j_1 + j_2$, with $J=|j_1 − j_2|, \codts , j_1 + j_2$.

### Clebsch-Gordan Coefficients

The CG coefficient $\ket j_1m_1;j_2m_2|JM\bra$ is computed via **Racah’s formula**. All quantum numbers are passed as doubled integers (e.g., $j=1/2 \rightarrow j_2=1$) to avoid floating‑point half‑integer issues.

```c
double clebsch_gordan(int j1_2, int m1_2, int j2_2, int m2_2,
                      int J_2, int M_2);
```

Returns 0 if selection rules are violated (triangle inequality, $m_1 + m_2 \neq M$, parity, etc.).

### Coupled State Construction

To build the coupled state $|j_1,j_2;J,M \bra$ in the uncoupled product basis {$|j_1,m_1 \bra \otimes |j_2,m_2 \bra}:

```c
cvector_t *couple_states(int j1_2, int j2_2, int J_2, int M_2);
```

The vector length is $(2j_1 + 1)(2j_2 + 1)$ (i.e., $(j1_2 + 1) * (j2_2 + 1)$). The basis ordering is row‑major: index $i1 * (j2_2+1) + i2$, where $i1$ corresponds to $m_1= −j_1 + i_1$ (step 1 in doubled units, i.e., $m_1$ increases by 1/2 each step). The same convention as used in `lz_matrix`.

To list all allowed total $J$ values:

```c
int J2_list[10];
int count = couple_allowed_J(j1_2, j2_2, J2_list);
// J2_list[0..count-1] contains allowed J values (doubled)
```

The spherical harmonics are the eigenfunctions of $L^2$ and $L_z$:

$$
Y_{lm}(\theta, \phi) = \sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}} P_l^m(\cos\theta) e^{i m \phi}
$$

> This page previously showed a `spherical_harmonic()` implementation, but > that function isn't declared in `angular.h` - it more likely lives in > `core/special/spherical_harmonics.c` alongside the other special functions > (Legendre, Laguerre, Hermite, Bessel). Worth pulling `core/special/special.h` in a later batch to document it in the right place rather than here.

## Spin Angular Momentum

Spin is an intrinsic angular momentum with $s=1/2$ for electrons:

$$
S_z|\uparrow\rangle = +\frac{\hbar}{2}|\uparrow\rangle, \quad
S_z|\downarrow\rangle = -\frac{\hbar}{2}|\downarrow\rangle
$$

### Implementation

The Pauli matrices are flat, row-major `extern` arrays (`sigma[i*2+j]`), not a 2D array or a wrapping struct:

```c
extern complex_t sigma_x[4];
extern complex_t sigma_y[4];
extern complex_t sigma_z[4];
```

Rather than exposing a generic "apply this 2×2 matrix" function, each Pauli operator has its own dedicated in-place application to a 2-component spinor (a plain `cvector_t` of length 2 - there's no separate `spinor_t` type):

```c
void spin_sigma_x(cvector_t *spinor);
void spin_sigma_y(cvector_t *spinor);
void spin_sigma_z(cvector_t *spinor);
```

```c
cvector_t *spin_up = cvector_alloc(2);
spin_up->data[0] = c_one();
spin_up->data[1] = c_zero();

spin_sigma_z(spin_up); // eigenstate: unchanged up to the +1 eigenvalue
```

For expectation values ($\langle S_z\rangle = \frac{\hbar}{2}\langle\psi|\sigma_z|\psi\rangle$, etc.), combine this with `cvector_dot` / `cvector_expect` from [Linear Algebra Core](../internals/linalg.md) - there's no dedicated `spin_expectation_z`-style helper in this header.

## Spin-Orbit Coupling

The spin-orbit interaction and its hydrogen-specific implementation in [Fine Structure](fine_structure.md) - `hydrogen_spin_orbit_energy`, `spin_orbit_ls_expect`, and the independent coupling-based cross-check `spin_orbit_ls_expect_from_coupling` (which is itself built from the Clebsch-Gordan machinery on this page).

## Coupled Angular Momentum

For two angular momenta $J_1$ and $J_2$, the total angular momentum is:

$$
\mathbf{J}=\mathbf{J}_1 + \mathbf{J}_2, \qquad j=|j_1-j_2|,\,|j_1-j_2|+1,\, \ldots,\, j_1+j_2
$$

The Clebsch-Gordan coefficients connect the uncoupled and coupled bases. Quantum numbers throughout this API are **doubled integers** (`j1_2 = 2*j1`, etc.) so half-integer spin stays representable as `int`:

```c
/* <j1 m1; j2 m2 | J M>, via Racah's formula. Condon-Shortley phase convention.
   Returns 0 if triangle inequality is violated, m1+m2 != M, |m| > j for any
   of the three, or j1+j2+J isn't an integer. */
double clebsch_gordan(int j1_2, int m1_2, int j2_2, int m2_2, int J_2, int M_2);

// |j1,j2; J,M> expressed in the uncoupled product basis {|j1,m1> x |j2,m2>}
cvector_t *couple_states(int j1_2, int j2_2, int J_2, int M_2);

/* Enumerate allowed total-J values (doubled):
   J_2 = |j1_2-j2_2|, |j1_2-j2_2|+2, ..., j1_2+j2_2 */
int couple_allowed_J(int j1_2, int j2_2, int *J2_out);
```

Example - coupling orbital $l=1$ with spin $s=\frac12$ (matching the header's own doc-comment example):

```c
/* spin-1/2, m=+1/2 -> j=1, m=+1 (i.e. j1_2=1, m1_2=1)
   orbital   l=1, m=0 -> j=2, m=0 (i.e. j2_2=2, m2_2=0)
   total     J=3/2    -> J_2=3 */
double cg = clebsch_gordan(1, 1, 2, 0, 3, 1);
```

`spin_orbit_ls_expect_from_coupling` uses to cross-check the closed-form spin-orbit expectation value in [Fine Structure](fine_structure.md).

## Example: Two Spin-1/2 Particles

```c
// Triplet M=+1 (J=1, M=1) -> |up up>
cvector_t *triplet_p1 = couple_states(1, 1, 2, 2);
// Singlet M=0 (J=0, M=0) -> (|up down> - |down up>)/sqrt2
cvector_t *singlet = couple_states(1, 1, 0, 0);
```

See `eg_15_angular_coupling.c` for a complete demo.

## Gaunt Coefficients

The Gaunt coefficient appears in multipole expansions of electron‑electron interactions:
$$G_{lm;l′m′}^k = \int Y_{lm}^{∗}(\Omega) Y_{kq}(\Omega) Y{l′m′}(\Omega)d\Omega, \quad q=m − m′$$

The function returns the **normalised** coefficient (without the $4\pi$ factor):

```c
double gaunt_coefficient(int l, int m, int k, int lp, int mp);
```

It uses the relation to Clebsch-Gordan coefficients and is validated against numerical quadrature in `test_angular_coupling.c`.

## See Also

- [Fine Structure](fine_structure.md) - uses spin‑orbit coupling and Clebsch-Gordan coefficients.
- [Identical Particles](identical.md) - Slater determinants and exchange symmetry.
- [Special Functions](../internals/special.md) - spherical harmonics implementation.
