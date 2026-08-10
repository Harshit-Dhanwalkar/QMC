# Harmonic Oscillator

The quantum harmonic oscillator is one of the most important exactly solvable models in quantum mechanics. It describes a particle in a parabolic potential $V(x) = \frac{1}{2}m\omega^2 x^2$.

## The Schr$\"{o}$dinger Equation

The time-independent Schr$\"{o}$dinger equation:

$$
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2} + \frac{1}{2}m\omega^2 x^2 \psi = E\psi
$$

The exact solutions are:

$$
\psi_n(x) = \frac{1}{\sqrt{2^n n!}} \left(\frac{m\omega}{\pi\hbar}\right)^{1/4}
H_n\!\left(\sqrt{\frac{m\omega}{\hbar}}x\right)
\exp\left(-\frac{m\omega x^2}{2\hbar}\right)
$$

Where $H_n$ are the Hermite polynomials.

The energy levels are equally spaced:

$$
E_n = \hbar\omega\left(n + \frac{1}{2}\right), \quad n = 0, 1, 2, \ldots
$$

Key properties:

- **Zero-point energy**: $E_0 = \frac{1}{2}\hbar\omega \neq 0$
- **Energy spacing**: $\Delta E = \hbar\omega$ (independent of n)
- **Parity**: Alternates between even ($n$ even) and odd ($n$ odd)
- **Number of nodes**: $n$ nodes

## Numerical Solution (Tridiagonal Diagonalization)

The finite-difference Hamiltonian for a 1D potential is real, symmetric, and tridiagonal - building a full dense `cmatrix_t` for it is wasted work. `core/ode/crank_nicolson.h` provides a shared builder for the `diag`/`offdiag` arrays, and `core/linalg/tridiag_eigh.h` diagonalizes them directly in $O(N^2)$ via the implicit QL algorithm (Wilkinson shift):

```c
void build_tridiagonal_hamiltonian(const double *x, const double *V, int N,
                                   double dx, double hbar_sq_2m, double *diag,
                                   double *offdiag);

eigen_t *tridiag_eigh(const double *diag, const double *offdiag, int n);
```

```c
double V[N];
for (int i = 0; i < N; i++)
    V[i] = 0.5 * m * omega * omega * x[i] * x[i];
```

> `physics/potentials.h` also provides `V_harmonic(x, params)` with `params: double *omega` - but note it has **no mass term**, computing `0.5*omega^2*x^2` directly (consistent with this project's natural-units convention, $\hbar=m=1$). If `m` isn't 1 in your setup, build the array manually as above rather than going through `V_harmonic` + `potential_array`. See [1D Potentials](potentials_1d.md).

```c
double diag[N], offdiag[N - 1];
build_tridiagonal_hamiltonian(x, V, N, dx, HBAR_SQ / (2.0 * m), diag, offdiag);

eigen_t *eig = tridiag_eigh(diag, offdiag, N);
```

`eig->eigenvalues` is ascending; `eig->eigenvectors` holds each eigenstate as a column of an $N \times N$ `cmatrix_t`. Extract an individual state with `core/utils.h`'s column helper:

```c
cvector_t *cvector_from_matrix_column(const cmatrix_t *m, int col);

cvector_t *psi_0 = cvector_from_matrix_column(eig->eigenvectors, 0);
```

This replaced an earlier dense-Jacobi approach (`cmatrix_eigh` on a full complex `cmatrix_t`) that ignored the matrix's tridiagonal structure - a meaningful correctness/performance change for every 1D potential module, not just the harmonic oscillator.

## Running the Example

```sh
./build/eg_02_harmonic
```

This generates:

- `harmonic_potential.dat` - the potential
- `harmonic_psi_0.dat` to `harmonic_psi_3.dat` - wavefunctions
- `harmonic_energies.dat` - energy levels
- `harmonic_plots.png` - visualization

### Verification

> Pending: needs actual numerical output compared against the analytic $E_n = \hbar\omega(n+1/2)$ spectrum to fill in properly - happy to add this once you share the test output or `test_tridiag.c` results.

The wavefunctions match the Hermite-Gauss form. The ground state is a Gaussian, while the first excited state is the odd-parity function $x\exp(−x^2/2)$.

### Expectation Values

For the harmonic oscillator:

$$
\langle x^2 \rangle =\frac{\hbar}{2m\omega}(2n+1),\, \langle p^2 \rangle=\frac{\hbar m \omega}{2}(2n+1)
$$

The uncertainty product saturates the Heisenberg bound:

$$
\Delta x \Delta p = \hbar \left(n+\frac{1}{2}\right)
$$

For the ground state, $\Delta x\Delta p=\hbar/2$, the minimum possible value.
