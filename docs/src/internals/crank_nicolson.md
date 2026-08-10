# Crank-Nicolson Solver

Time evolution of the 1D time-dependent Schrödinger equation for a tridiagonal Hamiltonian (finite-difference kinetic term plus a diagonal potential, optionally with a complex absorbing potential). Crank-Nicolson is unconditionally stable and unitary (to numerical precision) for real potentials, since it applies the Cayley form of the propagator:

$$
\left(I + \frac{i\,dt}{2}H\right)\psi_{n+1} = \left(I - \frac{i\,dt}{2}H\right)\psi_n
$$

Solved at each step via a complex tridiagonal Thomas algorithm - $O(N)$ per
step rather than $O(N^3)$ for a general linear solve.

Defined in `core/ode/crank_nicolson.h`.

## Building the Hamiltonian

```c
void build_tridiagonal_hamiltonian(const double *x, const double *V, int N,
                                   double dx, double hbar_sq_2m, double *diag,
                                   double *offdiag);
```

Discretizes $H = -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x)$ into `diag[i] = 2 * coeff + V[i]`, `offdiag[i] = -coeff`, where `coeff = hbar_sq_2m / dx^2`. This builder is shared with static diagonalization (see [Linear Algebra Core](linalg.md) - `tridiag_eigh` takes exactly this `diag`/`offdiag` pair), so the same discretization is used when finding eigenstates or propagating in time.

## Time-independent step

```c
int crank_nicolson_step(const double *diag, const double *offdiag, double dt,
                        cvector_t *psi);
```

One step of unitary evolution for a real, time-independent, diagonal-plus-offdiagonal Hamiltonian (no CAP).

## General step (complex diagonal / CAP)

```c
int crank_nicolson_step_general(const complex_t *diag, const double *offdiag,
                                double dt, cvector_t *psi);
```

Allows a complex diagonal, `diag[i] = V(x_i) [+ kinetic on-site term] - iΓ(x_i)`, where $\Gamma \geq 0$ is an optional complex absorbing potential (CAP). Pass an all-zero imaginary part to recover ordinary unitary evolution. The off-diagonal stays real, since CAP and any local potential only ever enter the Hamiltonian diagonally.

## Time-dependent potentials

```c
typedef double (*potential_time_fn)(double x, double t, void *params);

void build_tridiagonal_hamiltonian_time_dependent(
    const double *x, int N, double dx, double hbar_sq_2m, potential_time_fn V,
    void *params, double t, const double *absorb, complex_t *diag_out,
    double *offdiag_out);

int crank_nicolson_evolve_time_dependent(const double *x, int N, double dx,
                                         double hbar_sq_2m, potential_time_fn V,
                                         void *params, const double *absorb,
                                         cvector_t *psi, double t0, double dt,
                                         int steps);
```

`crank_nicolson_evolve_time_dependent` is the convenience driver: it steps `psi` forward `steps` times from `t0`, rebuilding $H$ at each step's **midpoint** time for second-order accuracy (rather than evaluating $V$ only at the step's start), and layering in `absorb` (maybe `NULL`) as a CAP.

## Complex absorbing potentials

```c
void cap_build_monomial(const double *x, int N, double width, double eta,
                        int power, double *absorb_out);
```

Riess-Meyer-style monomial CAP: $\Gamma$ ramps smoothly from 0 to `eta` over an absorbing layer of the given `width` at each edge of the grid ($\Gamma(x) = \eta \cdot ((x - x_{\text{layer start}})/\text{width})^{\text{power}}$), and is exactly 0 in the interior. `power` is typically 2 or 3. Used to prevent unphysical reflection off the grid boundary in scattering/tunneling simulations (see `eg_27_cap_tdse.c`).
