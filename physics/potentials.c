/*
V(x) library: infinite well, finite well, harmonic, barrier
*/

#include "potentials.h"
#include <math.h>

double V_infinite_well(double x, void *params) {
  double *a = (double *)params;
  if (fabs(x) < *a)
    return 0.0;
  return 1e10; /* Approximate infinity */
}

double V_finite_well(double x, void *params) {
  double *p = (double *)params;
  double a = p[0], V0 = p[1];
  if (fabs(x) < a)
    return -V0;
  return 0.0;
}

double V_harmonic(double x, void *params) {
  double *omega = (double *)params;
  return 0.5 * (*omega) * (*omega) * x * x;
}

double V_step(double x, void *params) {
  double *V0 = (double *)params;
  return (x >= 0) ? *V0 : 0.0;
}

double V_barrier(double x, void *params) {
  double *p = (double *)params;
  double a = p[0], V0 = p[1];
  if (x > 0 && x < a)
    return V0;
  return 0.0;
}

double V_coulomb(double r, void *params) {
  double *k = (double *)params; /* Coulomb constant */
  return -(*k) / r;
}

double V_yukawa(double r, void *params) {
  double *p = (double *)params;
  double g = p[0], mu = p[1];
  return -(g / r) * exp(-mu * r);
}

double V_morse(double x, void *params) {
  double *p = (double *)params;
  double D = p[0], a = p[1], x0 = p[2];
  double y = 1.0 - exp(-a * (x - x0));
  return D * y * y;
}

void potential_array(double *x, int n, potential_fn V, void *params,
                     double *V_out) {
  for (int i = 0; i < n; i++) {
    V_out[i] = V(x[i], params);
  }
}
