# Fine Structure Corrections

Fine structure corrections account for relativistic effects in single-electron systems, breaking the energy degeneracy with respect to total angular momentum $j$.

## Hamiltonian Corrections

The fine structure Hamiltonian consists of three terms:

$$H_{\text{fs}} = H_{\text{rel}} + H_{\text{SO}} + H_{\text{Darwin}}$$

### 1. Relativistic Kinetic Energy Correction

$$H_{\text{rel}} = -\frac{p^4}{8m^3c^2} \implies E_{\text{rel}}^{(1)} = -\frac{E_n^2}{2mc^2} \left( \frac{4n}{l + 1/2} - 3 \right)$$

### 2. Spin-Orbit Coupling

$$H_{\text{SO}} = \frac{e^2}{8\pi\varepsilon_0 m^2c^2 r^3} \mathbf{L} \cdot \mathbf{S}$$

### 3. Darwin Term ($l=0$)

$$H_{\text{Darwin}} = \frac{\pi\hbar^2 e^2}{2m^2c^2 (4\pi\varepsilon_0)} \delta^3(\mathbf{r})$$

## Combined Energy Shift

For hydrogenic atoms with atomic number $Z$:

$$E_{\text{fs}} = \frac{E_n \alpha^2 Z^2}{n} \left( \frac{1}{j + 1/2} - \frac{3}{4n} \right)$$

Where $\alpha \approx \frac{1}{137.036}$ is the fine-structure constant.

## Implementation

```c
#include "physics/fine_structure.h"

double fine_structure_shift(int n, int l, double j, int Z) {
    if (n <= 0 || l < 0 || l >= n || j < 0.5) return 0.0;

    double alpha = 1.0 / 137.035999139;
    double E_n = -13.605693 * Z * Z / (n * n); // eV

    double factor = (alpha * alpha * Z * Z) / n;
    double bracket = (1.0 / (j + 0.5)) - (3.0 / (4.0 * n));

    return E_n * factor * bracket;
}
```
