/*
Bessel functions (spherical harmonics)
*/
#include "special.h"
#include <math.h>
#include <stdlib.h>

/* Spherical Bessel j_l(x) using upward recurrence:
   j_{l+1} = (2l+1)/x * j_l - j_{l-1}
   with j_0 = sin(x)/x, j_1 = sin(x)/x^2 - cos(x)/x.
*/
double sph_bessel_j(int l, double x) {
  if (x < 1e-15) {
    if (l == 0)
      return 1.0;
    if (l == 1)
      return x / 3.0; // limit as x->0: j_l ~ x^l/(2l+1)!!
    return 0.0;
  }
  if (l == 0)
    return sin(x) / x;
  if (l == 1)
    return sin(x) / (x * x) - cos(x) / x;
  double j_prev2 = sin(x) / x;
  double j_prev1 = sin(x) / (x * x) - cos(x) / x;
  double j_cur = 0.0;
  for (int i = 1; i < l; i++) {
    j_cur = (2.0 * i + 1.0) / x * j_prev1 - j_prev2;
    j_prev2 = j_prev1;
    j_prev1 = j_cur;
  }
  return j_cur;
}

/* Spherical Neumann y_l(x) using recurrence with y_0 = -cos(x)/x, y_1 =
 * -cos(x)/x^2 - sin(x)/x */
double sph_bessel_y(int l, double x) {
  if (x < 1e-15)
    return -INFINITY; // singular
  if (l == 0)
    return -cos(x) / x;
  if (l == 1)
    return -cos(x) / (x * x) - sin(x) / x;
  double y_prev2 = -cos(x) / x;
  double y_prev1 = -cos(x) / (x * x) - sin(x) / x;
  double y_cur = 0.0;
  for (int i = 1; i < l; i++) {
    y_cur = (2.0 * i + 1.0) / x * y_prev1 - y_prev2;
    y_prev2 = y_prev1;
    y_prev1 = y_cur;
  }
  return y_cur;
}

/* Derivative: d/dx j_l(x) = (l/x) j_l(x) - j_{l+1}(x) */
double sph_bessel_j_deriv(int l, double x) {
  if (x < 1e-15) {
    if (l == 0)
      return 0.0;
    if (l == 1)
      return 1.0 / 3.0;
    return 0.0;
  }
  return (double)l / x * sph_bessel_j(l, x) - sph_bessel_j(l + 1, x);
}

void sph_bessel_array(int lmax, double x, double *j, double *y) {
  if (!j || !y || lmax < 0)
    return;
  if (x < 1e-15) {
    for (int l = 0; l <= lmax; l++) {
      j[l] = (l == 0) ? 1.0 : 0.0;
      y[l] = -INFINITY;
    }
    return;
  }
  j[0] = sin(x) / x;
  y[0] = -cos(x) / x;
  if (lmax >= 1) {
    j[1] = sin(x) / (x * x) - cos(x) / x;
    y[1] = -cos(x) / (x * x) - sin(x) / x;
  }
  for (int l = 2; l <= lmax; l++) {
    j[l] = (2.0 * l - 1.0) / x * j[l - 1] - j[l - 2];
    y[l] = (2.0 * l - 1.0) / x * y[l - 1] - y[l - 2];
  }
}
