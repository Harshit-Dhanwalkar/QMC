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

## Numerical Solution (Matrix Diagonalization)

The Hamiltonian is discretized on a grid using finite differences:

```c
// Build Hamiltonian matrix
double coeff = -hbar*hbar / (2*m*dx*dx);
for (int i = 0; i < N; i++) {
    double V = 0.5 * m * omega * omega * x[i] * x[i];
    CMAT(H, i, i) = c_real(-2.0*coeff + V);
    if (i > 0) CMAT(H, i, i-1) = c_real(coeff);
    if (i < N-1) CMAT(H, i, i+1) = c_real(coeff);
}

// Diagonalize
eigen_t *eig = cmatrix_eigh(H);
```

The eigenstates are obtained by diagonalizing the tridiagonal matrix.

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

// TODO

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
