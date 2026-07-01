# Numerov Integrator

The Numerov method is a fourth-order integrator for second-order ODEs of the form $y'' = f(x)y$ without first derivative terms - exactly the form of the 1D Schr$\"{o}$dinger equation.

## The Algorithm

For the Schr$\"{o}$dinger equation:

$$
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2} + V(x)\psi = E\psi
$$

rearrange as $\psi'' = f(x)\psi$ with $f(x) = \frac{2m}{\hbar^2}\left(V(x) - E\right)$.

Numerov's discretization (order $h^4$):

$$
\psi_{n+1} = \frac{2\psi_n(1 - 5h^2 f_n/12) - \psi_{n-1}(1 + h^2 f_{n-1}/12)}{1 + h^2 f_{n+1}/12}
$$

## Shooting Method

To find eigenvalues, we use the shooting method:

1. Guess energy $E$.
2. Integrate from left boundary ($\psi(0)=0$) to the right using Numerov.
3. The solution at the right boundary will be zero only if $E$ is an eigenvalue.
4. Use bisection to bracket and refine.

## Implementation

The core integrator is in `src/core/ode/numerov.c`:

```c
static double integrate_once(const numerov_params_t *p, double E,
                             cvector_t *psi) {
    int N = p->n;
    double h2 = p->dx * p->dx;
    double *f = malloc(N * sizeof *f);
    for (int i = 0; i < N; i++)
        f[i] = (p->V[i] - E) / p->hbar_sq_2m;

    psi->data[0].re = 0.0;
    psi->data[1].re = 1e-6;   // small initial slope
    for (int i = 1; i < N-1; i++) {
        double num = 2.0 * (1.0 - (5.0/12.0)*h2*f[i]) * psi->data[i].re
                   - (1.0 + (1.0/12.0)*h2*f[i-1]) * psi->data[i-1].re;
        double denom = 1.0 + (1.0/12.0)*h2*f[i+1];
        psi->data[i+1].re = num / denom;
    }
    free(f);
    return psi->data[N-1].re;   // value at right boundary
}
```

The bisection routine `numerov_shoot()` finds the energy where the boundary value crosses zero.

## Usage

```c
numerov_params_t p = { .x = x, .V = V, .n = N, .dx = dx, .hbar_sq_2m = 0.5 };
numerov_solution_t *sol = numerov_shoot(&p, 0.5, 1e-10);
printf("Energy: %f\n", sol->energy);
```

## Accuracy

Numerov is $O(h^4)$, so for typical grids (N \~ 1000) it gives $10^{−12}$ accuracy for eigenvalues.
