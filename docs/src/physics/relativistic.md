# Relativistic Quantum Mechanics

Relativistic quantum mechanics extends quantum theory to particles moving at speeds comparable to the speed of light.

## Klein-Gordon Equation

The relativistic wave equation for spin-0 particles:

$$
\left(\Box + \frac{m^2c^2}{\hbar^2}\right)\psi = 0
$$

Where $\Box = \frac{1}{c^2}\frac{\partial^2}{\partial t^2} - \nabla^2$, giving the energy-momentum relation $E^2 = p^2c^2 + m^2c^4$. This equation has both positive and negative energy solutions, leading to the concept of antiparticles.

`physics/relativistic.h` implements the 1D case, $\left(-\hbar^2c^2\frac{d^2}{dx^2} + m^2c^4\right)\psi = (E-V(x))^2\psi$, with **two different solvers** depending on whether the potential term needs to be treated exactly:

### Fast/bulk approximate spectrum

```c
eigen_t *klein_gordon_1d(double *x, int N, double *V, double m, double hbar,
                         double c);
```

This is the version that got switched from dense Jacobi diagonalization to `tridiag_eigh` (see [Linear Algebra Core](../internals/linalg.md)) as a correctness/performance fix. Because the true equation has $(E-V(x))^2$ on the right - which isn't a standard linear eigenvalue problem once $V(x)$ is spatially varying - this solver gives an _approximate_ bulk spectrum rather than an exact one; use `klein_gordon_1d_self_consistent` below when $V(x)$ varies enough that the approximation matters.

### Self-consistent solve (spatially-varying $V(x)$)

```c
typedef struct {
  double energy;
  cvector_t *psi;  // real-valued (im=0), dx-normalized
  int converged;   // 1 if fixed-point iteration converged within tol
  int iterations;
} klein_gordon_solution_t;

klein_gordon_solution_t *
klein_gordon_1d_self_consistent(double *x, int N, double *V, double m,
                                double hbar, double c, double E_guess,
                                double tol, int max_iter);

void klein_gordon_solution_free(klein_gordon_solution_t *sol);
```

Handles the $(E-V(x))^2$ nonlinearity properly via fixed-point iteration:

1. Guess $E$.
2. Build the tridiagonal eigenvalue problem $H(E)_{\text{diag}}[i] = 2\cdot\text{coeff} + m^2c^4 + 2E\,V(x_i) - V(x_i)^2$ and diagonalize.
3. Pick the eigenvalue $\lambda$ closest to the previous $E^2$.
4. $E_{\text{new}} = \text{sign}(E_{\text{guess}})\sqrt{\lambda}$; repeat until $|E_{\text{new}} - E| < \texttt{tol}$ or `max_iter` is hit.

`E_guess`'s **sign selects the energy branch** (positive-energy vs. negative-energy/antiparticle solutions) and its magnitude seeds which level gets tracked. Always check `sol->converged` - for strongly-varying $V(x)$ the iteration can exhaust `max_iter` without converging, and `sol->energy`/`psi` are still the last-iterate values regardless of whether it actually converged.

```c
klein_gordon_solution_t *sol = klein_gordon_1d_self_consistent(
    x, N, V, m, hbar, c, /*E_guess=*/m * c * c, 1e-10, 200);
if (!sol->converged) {
    // last-iterate values still populated
}
klein_gordon_solution_free(sol);
```

## Dirac Equation

The general relativistic wave equation for spin-1/2 particles in 3D is $i\hbar\,\partial_t\psi = (c\,\boldsymbol{\alpha}\cdot\mathbf{p} + \beta * mc^2)\psi$ with 4×4 $\boldsymbol{\alpha}$, $\beta$ matrices and a four-component spinor.

**What's actually implemented is the 1D reduction**, a 2-component equation using Pauli matrices instead of the full 4×4 Dirac matrices:

```c
// [ c * \sigma_z * p + m * c^2 * \sigma_x + V(x) ] * \phi = E * \phi
eigen_t *dirac_1d(double *x, int N, double *V, double m, double hbar, double c);
```

Returns an `eigen_t` whose eigenvectors are 2-component spinors. This uses `sigma_x`/`sigma_z` from [Angular Momentum & Spin](angular_momentum.md), not a separate 4×4 `gamma_matrix_t` - the previous version of this page showed a 4-component `dirac_spinor_t` with an analytic plane-wave constructor (`dirac_plane_wave`), but that isn't what `relativistic.h` declares. If a plane-wave/4-spinor path exists it's elsewhere and needs its own header before documenting.

## Non-Relativistic Limit

In the limit $v \ll c$, the Dirac equation reduces to the Schrödinger equation with spin-orbit coupling - see [Fine Structure](fine_structure.md) for the implemented spin-orbit terms.

### Applications

> The closed-form Dirac fine-structure formula ($E_{nj} = mc^2\left[1+\left(\frac{\alpha}{n-(j+1/2)+\sqrt{(j+1/2)^2-\alpha^2}}\right)^2\right]^{-1/2}$) isn't backed by a function in this header. It's kept here as physics background/a target for validating `dirac_1d`'s hydrogen-like spectrum against, not as documented API.
