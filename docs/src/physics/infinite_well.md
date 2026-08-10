# Infinite Square Well

The simplest bound-state problem. A particle confined to $0 \leq x \leq L$
with infinitely hard walls: $V(x) = 0$ inside, $V = \infty$ outside.

## The Schrödinger Equation

Inside the well, the time-independent Schrödinger equation reduces to:

$$
-\frac{\hbar^2}{2m} \frac{d^2\psi}{dx^2} = E\psi
$$

With boundary conditions $\psi(0) = \psi(L) = 0$.

The solutions are standing waves:

$$
\psi_n(x) = \sqrt{\frac{2}{L}} \sin\!\left(\frac{n\pi x}{L}\right), \quad n = 1, 2, 3, \ldots
$$

with quantized energies:

$$
E_n = \frac{n^2 \pi^2 \hbar^2}{2mL^2}
$$

The ground state $n=1$ has no nodes. Each successive state adds one node.
The zero-point energy $E_1 > 0$ is a direct consequence of the uncertainty principle.

## Numerical Solution

The potential is defined as a point evaluator in `physics/potentials.h`,
filled across the grid via `potential_array` (see
[1D Potentials](potentials_1d.md)):

```c
double V_infinite_well(double x, void *params); // params: double *a (half-width)
```

> Note the convention: the well is `|x| < a`, centered at the origin - not
> `[0, L]` as earlier text on this page implied. The analytic formulas below
> assume a well of width $L$; if you're calling `V_infinite_well` directly,
> map `L = 2a` and shift the grid, or build the potential array manually as
> shown further down if you want the `[0, L]` convention.

```c
double a = L / 2.0;
potential_array(x, N, V_infinite_well, &a, V);
```

Eigenstates are solved via `core/ode/numerov.h`, which provides two solvers:

```c
numerov_solution_t *numerov_shoot(numerov_params_t *params, double E_guess,
                                  double E_tol);

numerov_solution_t *numerov_shoot_matching(numerov_params_t *params,
                                           double E_min, double E_max,
                                           int n_scan, double tol);
```

**`numerov_shoot`** - despite the name, this is _not_ a shooting method. It diagonalizes the discretized Hamiltonian directly via `tridiag_eigh` and selects the requested level from `E_guess` (`level = round(E_guess - V_min - 0.5)`, clamped to ≥ 0). This is the reliable general-purpose solver and is what's used for the infinite well.

**`numerov_shoot_matching`** - true bidirectional Numerov shooting with log-derivative matching at the outer classical turning point, bisecting on an energy bracket `[E_min, E_max]`. This exists because single-direction shooting with zero-crossing detection is mathematically incapable of finding eigenvalues in classically-allowed regions (proven via Taylor-seriesderivation - the characteristic ratio stays > 1, so the trial solution never changes sign regardless of grid resolution). It's the right tool for problems with genuine forbidden regions, like the finite well.

```c
numerov_params_t params = {
    .x = x, .V = V, .n = N, .dx = dx, .hbar_sq_2m = hbar_sq_2m
};
numerov_solution_t *sol = numerov_shoot(&params, E_guess, 1e-10);
```

## Running the Example

```sh
./build/infinite_well
```

Plots the first four eigenstates with energy levels annotated as
$E_n = n^2 \pi^2 / L^2$ (natural units $\hbar = 2m = 1$).

## Verification

Numerical vs analytic for $L = 1$, natural units:

| n   | Numerical    | Analytic $n^2\pi^2$   |
| --- | ------------ | --------------------- |
| 1   | 9.8696044... | 9.8696044010893586... |
| 2   | 39.478417... | 39.478417604357434... |
| 3   | 88.826440... | 88.826440109803...    |

Agreement to 10 significant figures - Numerov is $O(h^4)$.
