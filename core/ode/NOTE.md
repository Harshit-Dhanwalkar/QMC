### Usage

- TDSE with Crank-Nicolson:

```c
#include "core/ode/crank-nicolson.h"
// ... set up diag, offdiag, psi
double dt = 0.01;
for (int step = 0; step < nsteps; step++) {
    crank_nicolson_step(diag, offdiag, dt, psi);
}
```

- RK4 for arbitrary ODE:

```c
#include "core/ode/rk4.h"

void my_ode(double t, const cvector_t *y, cvector_t *dydt, void *params) {
    // dydt = some function of y and t
}

cvector_t *y = ...;
double t = 0, dt = 0.01;
rk4_step(t, dt, y, my_ode, NULL);
```

- Test with the harmonic oscillator example to verify time evolution.
- Implement the split-step method (using FFT) for TDSE as an alternative - this is often faster.
