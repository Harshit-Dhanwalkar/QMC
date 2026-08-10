# Time-Dependent Schrödinger Equation Integrators

(TISE/TDSE Convenience Layer)

Provides real- and imaginary-time wave function evolution solvers.

## Time Evolution Operator

$$
\psi(x, t + \Delta t) = e^{-\frac{i H \Delta t}{\hbar}} \psi(x, t)
$$

## Integration Algorithms

### 1. Crank-Nicolson Scheme (Unitary & Implicit)

$$
\left( I + \frac{i \Delta t}{2\hbar} H \right) \psi(x, t + \Delta t) = \left( I - \frac{i \Delta t}{2\hbar} H \right) \psi(x, t)
$$

Solves tridiagonal matrix system using Thomas algorithm.

### 2. Imaginary Time Propagation (Ground State Calculation)

Substituting $t \to -i\tau$:

$$
\psi(x, \tau + \Delta \tau) = e^{-\frac{H \Delta \tau}{\hbar}} \psi(x, \tau)
$$

Normalizing $\psi$ at each step projects out excited states, converging to the ground state energy.

## Time-independent (TISE)

```c
/* Solve 1D TISE by matrix diagonalization (finite difference).
   Eigenvectors are stored as columns of the eigenvectors matrix.
   */
eigen_t *solve_tise_matrix(double *x, int n, double dx, double hbar_sq_2m,
                           potential_fn V, void *params);
```

NOTE: eigenvectors are not normalized to `dx` - call `wavefunction_normalize()` (see [Wave functions](wavefunctions.md)) after extracting one.

Takes a [1D-potentials-style](potentials_1d.md) point evaluator directly (consistent with `central_potential_radial_solve` and the `wkb.h`/ `scattering.h` functions), builds the tridiagonal Hamiltonian internally, and diagonalizes via `tridiag_eigh` - the same pattern used by-hand in [Harmonic Oscillator](harmonic_oscillator.md), just wrapped into one call.

```c
// Shooting method (Numerov) for one specific bound state
numerov_solution_t *solve_tise_shoot(numerov_params_t *params, double E_guess,
                                     double E_tol);
```

Thin wrapper - check whether this forwards to `numerov_shoot` or `numerov_shoot_matching` (see [Numerov Integrator](../internals/numerov.md); the two behave very differently) before relying on it for a potential with genuine classically-forbidden regions.

## Time-dependent (TDSE)

```c
/* Crank-Nicolson evolution. H tridiagonal (diag, offdiag).
    Returns 0 on success. */
int evolve_tdse_crank(double *diag, double *offdiag, int n, cvector_t *psi,
                      double dt, int steps);
```

A simplified front door onto [`crank_nicolson_step`](../internals/crank_nicolson.md) repeated `steps` times - note this signature takes real `diag`/`offdiag` (no CAP, no time-dependence), so it's the plain unitary case only; for CAP or time-dependent potentials, go directly to `crank_nicolson_evolve_time_dependent`.

```c
/* Split-step Fourier method: propagates via FFT to momentum space and back.
   Only for potentials that are functions of x alone. */
int evolve_tdse_split_step(cvector_t *psi, double *x, double *V, int n,
                           double dx, double dt, int steps, double hbar,
                           double mass);
```

The 1D counterpart to [SOFT](soft.md)'s 2D/3D split-operator propagation - same Strang-splitting idea, one dimension.

## TODO
When to use this vs. lower-level modules directly

See `schrodinger.h` for a quick solve where the defaults are fine; reach for the internals pages directly when need CAP, time-dependent potentials, the shooting-vs-diagonalization choice made explicit, or control over normalization timing.
