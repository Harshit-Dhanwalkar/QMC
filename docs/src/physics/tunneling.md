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

For $E < V_0$, the wavefunction decays exponentially in the barrier ($\kappa = \sqrt{2m(V_0-E)}/\hbar$) and $R = 1$ (complete reflection, with exponential penetration).

## Tunneling Barrier

For a rectangular barrier of height $V_0$ and width $a$, the **exact** transmission probability (from matching wavefunctions and derivatives at both edges, not an approximation) is:

For $E < V_0$:

$$
T = \frac{1}{1 + \frac{V_0^2}{4E(V_0 - E)} \sinh^2(\kappa a)}, \quad \kappa = \sqrt{2m(V_0-E)}/\hbar
$$

For $E > V_0$:

$$
T = \frac{1}{1 + \frac{V_0^2}{4E(E - V_0)} \sin^2(ka)}, \quad k = \sqrt{2m(E-V_0)}/\hbar
$$

> **Implementation status:** this page previously showed a `barrier_transmission()` C function implementing the exact formulas above, but no header reviewed so far (`potentials.h`, `wkb.h`, `scattering.h`) declares a matching function - it may live in a module not yet reviewed, or the exact closed-form result may currently only exist as physics reference here rather than as callable code. What **is** implemented and confirmed is the WKB _approximation_ to this same barrier, in [WKB Approximation](wkb.md):
>
> ```c
> /* T = exp(-2*kappa*width), kappa = sqrt((V0-E)/hbar_sq_2m).
>    Returns 1.0 if V0 <= E. */
> double wkb_tunneling_rectangular(double E, double V0, double width,
>                                  double hbar_sq_2m);
> ```
>
> This is the leading-exponential approximation to the exact $\sinh^2$ formula above - good when $\kappa a \gg 1$ (thick/high barrier), where $\sinh^2(\kappa a) \approx \frac14 e^{2\kappa a}$ dominates. The two should agree closely in that regime and diverge for thin/low barriers, which would be a good numerical cross-check once the exact version's actual function (if any) is located.

## Resonant Tunneling

When $E>V_0$, the transmission shows resonant peaks when $\sin(ka)=0$, i.e. when $ka=n\pi$. This corresponds to constructive interference inside the barrier.
