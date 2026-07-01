# Infinite Square Well

The simplest bound-state problem. A particle confined to $0 \leq x \leq L$
with infinitely hard walls: $V(x) = 0$ inside, $V = \infty$ outside.

## The Schrödinger Equation

Inside the well, the time-independent Schrödinger equation reduces to:

$$
-\frac{\hbar^2}{2m} \frac{d^2\psi}{dx^2} = E\psi
$$

with boundary conditions $\psi(0) = \psi(L) = 0$.

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

We solve this via the **Numerov method** on a discrete grid. See
[Numerov Integrator](../internals/numerov.md) for the algorithm detail.

The potential is defined in `src/physics/potentials.c`:

```c
void potential_infinite_well(double *V, const double *x, int N, double L) {
    for (int i = 0; i < N; i++)
        V[i] = (x[i] >= 0.0 && x[i] <= L) ? 0.0 : 1e30;
}
```

Eigenvalues are found by bisection in `src/core/numerov.c` via `shoot_eigenvalue()`.
The solver brackets each level between consecutive analytic values, then bisects
until `|E_hi - E_lo| < 1e-10`.

## Running the Example

```sh
./build/infinite_well
```

Plots the first four eigenstates with energy levels annotated as
$E_n = n^2 \pi^2 / L^2$ (natural units $\hbar = 2m = 1$).

## Verification

Numerical vs analytic for $L = 1$, natural units:

| n   | Numerical    | Analytic $n^2\pi^2$ |
| --- | ------------ | ----------------------- |
| 1   | 9.8696044... | 9.8696044010893586...   |
| 2   | 39.478417... | 39.478417604357434...   |
| 3   | 88.826440... | 88.826440109803...      |

Agreement to 10 significant figures - Numerov is $O(h^4)$.
