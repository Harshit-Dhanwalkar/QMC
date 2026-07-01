# Step Potential & Tunneling

A particle incident on a potential step or barrier demonstrates quantum tunneling - a purely quantum effect where particles can pass through classically forbidden regions.

## Step Potential

Consider a potential step:

$$
V(x) = \begin{cases}
0 & x < 0 \\
V_0 & x \geq 0
\end{cases}
$$

For $E > V_0$, the solutions are:

**Region I** ($x < 0$): $\psi(x) = e^{ikx} + re^{-ikx}$, $k = \sqrt{2mE}/\hbar$

**Region II** ($x \geq 0$): $\psi(x) = te^{iqx}$, $q = \sqrt{2m(E - V_0)}/\hbar$

The reflection and transmission coefficients:

$$
R = \left|\frac{k - q}{k + q}\right|^2, \quad T = \frac{4kq}{(k + q)^2}
$$

For $E < V_0$, the wavefunction decays exponentially in the barrier:

**Region II**: $\psi(x) = te^{-\kappa x}$, $\kappa = \sqrt{2m(V_0 - E)}/\hbar$

The reflection coefficient becomes $R = 1$ (complete reflection), but with exponential penetration.

## Tunneling Barrier

For a rectangular barrier of height $V_0$ and width $a$:

$$
V(x) = \begin{cases}
0 & x < 0 \\
V_0 & 0 \leq x \leq a \\
0 & x > a
\end{cases}
$$

The transmission probability for $E < V_0$ is:

$$
T = \frac{1}{1 + \frac{V_0^2}{4E(V_0 - E)} \sinh^2(\kappa a)}
$$

where $\kappa = \sqrt{2m(V_0 - E)}/\hbar$.

For $E > V_0$:

$$
T = \frac{1}{1 + \frac{V_0^2}{4E(E - V_0)} \sin^2(ka)}
$$

with $k = \sqrt{2m(E - V_0)}/\hbar$.

## Implementation

```c
double barrier_transmission(double E, double V0, double a,
                           double m, double hbar) {
    if (E < 0) return 0.0;
    double kappa = sqrt(2.0 * m * (V0 - E)) / hbar;
    if (E < V0) {
        double sinh_term = sinh(kappa * a);
        double denom = 1.0 + (V0*V0) / (4.0*E*(V0-E)) * sinh_term*sinh_term;
        return 1.0 / denom;
    } else {
        double k = sqrt(2.0 * m * (E - V0)) / hbar;
        double sin_term = sin(k * a);
        double denom = 1.0 + (V0*V0) / (4.0*E*(E-V0)) * sin_term*sin_term;
        return 1.0 / denom;
    }
}
```

## Resonant Tunneling

When $E>V_0$, the transmission shows resonant peaks when $\sin(ka)=0$, i.e., when $ka=n\pi k$. This corresponds to constructive interference inside the barrier.
