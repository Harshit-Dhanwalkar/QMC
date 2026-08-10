# Finite Square Well

A particle trapped in a finite potential well: $V(x) = -V_0$ for $|x| < L/2$ and $V(x) = 0$ outside. Unlike the infinite well, the wavefunction can tunnel into the classically forbidden regions.

## The Schr$\ddot{o}$dinger Equation

The finite square well potential:

$$
V(x) = \begin{cases}
-V_0 & |x| < L/2 \\
0 & |x| \geq L/2
\end{cases}
$$

The time-independent Schrödinger equation must be solved separately in three regions:

**Region I** ($x < -L/2$): $\psi(x) = Ae^{\kappa x}$

**Region II** ($-L/2 \leq x \leq L/2$): $\psi(x) = B\cos(kx) + C\sin(kx)$

**Region III** ($x > L/2$): $\psi(x) = De^{-\kappa x}$

Where:

$$
k = \frac{\sqrt{2m(E + V_0)}}{\hbar}, \quad \kappa = \frac{\sqrt{-2mE}}{\hbar}
$$

The quantized energy levels are found by matching $\psi$ and $d\psi/dx$ at the boundaries, yielding transcendental equations:

**Even parity**: $\kappa = k\tan(kL/2)$

**Odd parity**: $\kappa = -k\cot(kL/2)$

The number of bound states is finite, given by:

$$
N = \left\lfloor \frac{L}{\pi\hbar}\sqrt{2mV_0} \right\rfloor + 1
$$

## Numerical Solution

The potential is defined as a point evaluator (see
[1D Potentials](potentials_1d.md)):

```c
double V_finite_well(double x, void *params); // params: struct { double a, V0; }
```

Here `a` is the **half-width** (well spans `|x| < a`, so `a = L/2` in the
notation above):

```c
struct { double a, V0; } params = { .a = L / 2.0, .V0 = V0 };
double *V = malloc(N * sizeof(double));
potential_array(x, N, V_finite_well, &params, V);
```

Because the wavefunction has non-trivial amplitude in the classically forbidden regions outside the well, this is exactly the case `numerov_shoot_matching` (bidirectional shooting with log-derivative matching, see [Numerov Integrator](../internals/numerov.md)) exists for - the single-direction / matrix-diagonalization approach used for the infinite well isn't the right tool here.

## Running the Example

```sh
./build/finite_well
```
