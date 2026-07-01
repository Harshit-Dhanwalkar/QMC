Usage in examples

```c 
#include "../core/fft/fft.h"
// ...
cvector_t *psi_x = ...; // your wavefunction on a grid
double dx = ...;
cvector_t *psi_k = position_to_momentum(psi_x, dx);

// Now psi_k->data contains the momentum-space wavefunction.
// The corresponding k values can be generated as:
int N = psi_k->n;
double dk = 2.0 * M_PI / (N * dx);
for (int i = 0; i < N; i++) {
    double k = (i < N/2) ? i*dk : (i - N)*dk;   // after fft_shift
}
// Then you can compute expectation values, plot, etc.

cvector_free(psi_k);
```
