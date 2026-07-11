#ifndef QMC_RK4_H
#define QMC_RK4_H

#include "../vector.h"

/* Runge-Kutta 4th order integrator for complex ODEs.
   System: d\phi/dt = f(\phi, t).
*/

typedef void (*ode_func_t)(double t, const cvector_t *y, cvector_t *dydt,
                           void *params);

/* Perform one RK4 step from time t to t+dt.
   y: input/output vector (complex)
   params: additional parameters for f.
   Returns 0 on success.
*/
int rk4_step(double t, double dt, cvector_t *y, ode_func_t f, void *params);

#endif
