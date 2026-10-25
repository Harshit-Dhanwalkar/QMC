/*
 * 1D and radial potential library
 *
 * Provides potential callbacks used by 1D Schrodinger solvers and radial atomic
 * solvers
 *
 * Included potentials:
 *   V_infinite_well(x, &a)       0 inside |x| < a, ~\infty outside
 *   V_finite_well(x, {a, V0})    -V0 inside |x| < a, 0 outside
 *   V_harmonic(x, &omega)        (1/2) * \omega^2 * x^2
 *   V_step(x, &V0)               V0 for x >= 0, 0 otherwise
 *   V_barrier(x, {a, V0})        V0 for 0 < x < a, 0 otherwise
 *   V_coulomb(r, &k)             -k / r  (radial)
 *   V_yukawa(r, {g, mu})        -(g/r) * \exp(-\mu * r)  (radial)
 *   V_morse(x, {D, a, x0})       D * (1 - \exp(-a * (x-x0)))^2
 *
 * potential_array() evaluates any one of these on a grid: it is a thin
 * loop calling V(x[i], params) for i in [0, n), writing V_out[i]
 */

#include "potentials.h"
#include <math.h>

double V_infinite_well(double x, const void *params) {
  const double *a = (const double *)params;
  if (fabs(x) < *a) {
    return 0.0;
  }

  return 1e10; // Approximate infinity
}

double V_finite_well(double x, const void *params) {
  const double *p = (const double *)params;
  double a = p[0];
  if (fabs(x) < a) {
    return -p[1];
  }

  return 0.0;
}

double V_harmonic(double x, const void *params) {
  const double *omega = (const double *)params;

  return 0.5 * (*omega) * (*omega) * x * x;
}

double V_step(double x, const void *params) {
  const double *V0 = (const double *)params;

  return (x >= 0) ? *V0 : 0.0;
}

double V_barrier(double x, const void *params) {
  const double *p = (const double *)params;
  double a = p[0];
  if (x > 0 && x < a) {
    return p[1];
  }

  return 0.0;
}

double V_coulomb(double r, const void *params) {
  const double *k = (const double *)params; // Coulomb constant

  return -(*k) / r;
}

double V_yukawa(double r, const void *params) {
  const double *p = (const double *)params;
  double g = p[0];
  double mu = p[1];

  return -(g / r) * exp(-mu * r);
}

double V_morse(double x, const void *params) {
  const double *p = (const double *)params;
  double D = p[0];
  double a = p[1];
  double x0 = p[2];
  double y = 1.0 - exp(-a * (x - x0));

  return D * y * y;
}

void potential_array(double *x, int n, potential_fn V, const void *params,
                     double *V_out) {
  for (int i = 0; i < n; i++) {
    V_out[i] = V(x[i], params);
  }
}
