# Finite Square Well

A particle trapped in a finite potential well: $V(x) = -V_0$ for $|x| < L/2$ and $V(x) = 0$ outside. Unlike the infinite well, the wavefunction can tunnel into the classically forbidden regions.

## The Schrödinger Equation

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

where:

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

The finite well is implemented in `src/physics/potentials.c`:

```c
void potential_finite_well(double *V, const double *x, int N,
                          double V0, double L) {
    for (int i = 0; i < N; i++) {
        V[i] = (fabs(x[i]) < L/2.0) ? -V0 : 0.0;
    }
}
```

## Running the Example

```sh
./build/finite_well
```
