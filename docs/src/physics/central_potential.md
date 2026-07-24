# Central Potentials

Central potentials are spherically symmetric potentials $V(\mathbf{r}) = V(r)$ that depend only on the radial distance from the origin.

## Radial Schrödinger Equation

In three dimensions, the wave function separates into radial and angular parts:

$$
\psi(r, \theta, \phi) = R_{nl}(r) Y_{lm}(\theta, \phi) = \frac{u_{nl}(r)}{r} Y_{lm}(\theta, \phi)
$$

The reduced radial wave function $u(r) = r R(r)$ satisfies the 1D effective Schrödinger equation:

$$
-\frac{\hbar^2}{2m} \frac{d^2 u}{dr^2} + V_{\text{eff}}(r) u(r) = E u(r)
$$

Where the effective potential includes the centrifugal barrier:

$$
V_{\text{eff}}(r) = V(r) + \frac{\hbar^2 l(l+1)}{2m r^2}
$$

## Implementation

Radial equations are solved on a discrete grid starting near $r \to 0$ with boundary condition $u(0) = 0$.

```c
#include "physics/central_potential.h"

// Computes effective central potential on a grid
void central_effective_potential(double *V_eff, const double *V, const double *r,
                                 int N, int l, double m, double hbar) {
    for (int i = 0; i < N; i++) {
        if (r[i] <= 0.0) {
            V_eff[i] = V[i];
            continue;
        }
        double centrifugal = (hbar * hbar * l * (l + 1)) / (2.0 * m * r[i] * r[i]);
        V_eff[i] = V[i] + centrifugal;
    }
}
```
