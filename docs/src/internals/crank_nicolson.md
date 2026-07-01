# Crank-Nicolson Solver

The Crank-Nicolson method is an unconditionally stable and unitary time propagator for the time-dependent Schrödinger equation (TDSE).

## Method

The TDSE: $i\hbar \frac{\partial\psi}{\partial t} = H\psi$.

Discretize time with step $\Delta t$:

$$
\frac{\psi^{n+1} - \psi^n}{\Delta t} = -\frac{i}{\hbar} H \frac{\psi^{n+1} + \psi^n}{2}
$$

Rearranging gives:

$$
\left(I + \frac{i\Delta t}{2\hbar}H\right)\psi^{n+1} = \left(I - \frac{i\Delta t}{2\hbar}H\right)\psi^n
$$

For a 1D system with tridiagonal Hamiltonian, this becomes a complex tridiagonal system solvable by the Thomas algorithm.

## Implementation

The core function is `crank_nicolson_step()` in `src/core/ode/crank_nicolson.c`:

```c
int crank_nicolson_step(const double *diag, const double *offdiag, double dt,
                        cvector_t *psi) {
    int N = psi->n;
    // Build A = I + i*dt/2*H (tridiagonal)
    complex_t *a = malloc(...); // lower
    complex_t *b = malloc(...); // diagonal
    complex_t *c = malloc(...); // upper
    // Compute RHS = (I - i*dt/2*H) * psi
    // Solve A * psi_new = RHS using Thomas algorithm
    // Overwrite psi with solution
}
```

The Thomas algorithm for complex tridiagonal systems is implemented internally.

## Building the Hamiltonian

The function `build_tridiagonal_hamiltonian` constructs the diagonal and off-diagonal arrays:

```c
void build_tridiagonal_hamiltonian(const double *x, const double *V, int N,
                                   double dx, double hbar_sq_2m,
                                   double *diag, double *offdiag) {
    double coeff = hbar_sq_2m / (dx*dx);
    for (int i = 0; i < N; i++) {
        diag[i] = 2.0*coeff + V[i];
        if (i < N-1) offdiag[i] = -coeff;
    }
}
```

### Usage Example

```c
double *x = linspace(0, L, N);
double *V = calloc(N, sizeof(double));
double *diag = malloc(N * sizeof(double));
double *offdiag = malloc((N-1) * sizeof(double));
build_tridiagonal_hamiltonian(x, V, N, dx, HBAR_2M, diag, offdiag);

cvector_t *psi = initial_wavepacket(...);
for (int step = 0; step < n_steps; step++) {
    crank_nicolson_step(diag, offdiag, dt, psi);
    // save or plot psi at desired times
}
```

## Properties

- Unitary: preserves norm (probability) exactly.
- Unconditionally stable: no restriction on $\Delta t$ (though accuracy improves with smaller steps).
- Second-order accurate in time and space (when used with second-order spatial discretization).
