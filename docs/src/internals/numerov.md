# Numerov Integrator

The Numerov method solves the 1D time-independent Schr$\ddot{o}$dinger equation

$$
-\frac{\hbar^2}{2m}\frac{d^2\phi}{dx^2} + V(x)\phi = E\phi
$$

By exploiting the absence of a first-derivative term: for an ODE of the form $\phi''(x) = f(x)\phi(x)$, Numerov gives $O(h^4)$ local accuracy from a three-point recurrence, better than the $O(h^2)$ of a generic finite-difference scheme at the same cost.

Defined in `core/ode/numerov.h`.

## Grid & Parameters

```c
typedef struct {
  double *x;         // Position grid
  double *V;         // Potential array V(x)
  int n;             // Grid points
  double dx;         // Grid spacing
  double hbar_sq_2m; // \hbar^2/(2m) in problem units
} numerov_params_t;

typedef struct {
  double energy;
  cvector_t *psi;
} numerov_solution_t;
```

## Two solvers, two purposes

QMC exposes two entry points that both produce a `numerov_solution_t`, but they work in fundamentally different ways:

### `numerov_shoot` - matrix diagonalization

```c
numerov_solution_t *numerov_shoot(numerov_params_t *params, double E_guess,
                                  double E_tol);
```

`E_guess` is used only to pick a target level index (`level = round(E_guess - V_min - 0.5)`, clamped to ≥ 0); the eigenvalue and eigenvector are then obtained by discretizing the Hamiltonian into `diag`/`offdiag` arrays and calling `tridiag_eigh` (see [Linear Algebra Core](linalg.md)). This is the default, reliable path for problems without classically forbidden regions to match across - e.g. the infinite square well.

### `numerov_shoot_matching` - true shooting

```c
numerov_solution_t *numerov_shoot_matching(numerov_params_t *params,
                                           double E_min, double E_max,
                                           int n_scan, double tol);
```

Bidirectional Numerov integration with log-derivative matching at the outer classical turning point, bracketed on `[E_min, E_max]` and refined by bisection (`n_scan` points used to find a sign change before bisecting). Validated to roughly $10^{-9}$–$10^{-11}$ against the harmonic oscillator's known spectrum.

This exists because naive single-direction shooting - integrate outward from one boundary and look for a zero-crossing at some matching point - cannot detect eigenvalues in a classically-allowed region: a Taylor-series analysis of the recurrence shows the characteristic ratio stays above 1, so the trial solution never changes sign there regardless of step size or grid resolution. Bidirectional matching (integrate from both ends and match log-derivatives at a turning point) sidesteps the problem entirely.

## Raw integration

```c
void numerov_integrate(const numerov_params_t *params, double E,
                       cvector_t *psi);
```

Forward Numerov integration at a fixed energy, seeded with $\psi_0 = 0$ (Dirichlet), $\psi_1 = 10^{-8}$. Used internally by both solvers above; also useful standalone for plotting a trial wavefunction at an arbitrary (possibly non-eigen) energy.

## Cleanup

```c
void numerov_solution_free(numerov_solution_t *sol);
```
