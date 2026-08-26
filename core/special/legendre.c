/*
Legendre polynomials (angular momentum)
*/
#include "special.h"
#include <math.h>
#include <stdlib.h>

/* Legendre polynomials P_l(x) using recurrence:
 *  (l+1)P_{l+1} = (2l+1)x P_l - l P_{l-1}
 * P_0=1, P_1=x
 */
double legendre(int l, double x) {
  if (l < 0) {
    return 0.0;
  }
  if (l == 0) {
    return 1.0;
  }
  if (l == 1) {
    return x;
  }

  double P_prev2 = 1.0;
  double P_prev1 = x;
  double P_cur = 0.0;
  for (int i = 1; i < l; i++) {
    P_cur = ((2.0 * i + 1.0) * x * P_prev1 - i * P_prev2) / (i + 1.0);
    P_prev2 = P_prev1;
    P_prev1 = P_cur;
  }

  return P_cur;
}

/* Associated Legendre P_l^m(x) using recurrence or relation.
 * Use recurrence from Numerical Recipes:
 *   P_l^m = (-1)^m * (1 - x^2)^{m/2} d^m/dx^m P_l(x)
 * Implement using recursion in l for fixed m.
   HACK: for small l, compute via derivative of P_l.
*/
double assoc_legendre(int l, int m, double x) {
  if (m < 0 || m > l) {
    return 0.0;
  }
  if (l == 0 && m == 0) {
    return 1.0;
  }

  // Use recurrence on l for fixed m.
  // Base: P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}
  // P_{m+1}^m = (2m+1) x P_m^m
  // General: (l-m) P_l^m = (2l-1) x P_{l-1}^m - (l+m-1) P_{l-2}^m
  if (m == 0) {
    return legendre(l, x);
  }

  // Compute P_m^m
  double P_mm = 1.0;
  for (int i = 1; i <= m; i++) {
    P_mm *= -(2.0 * i - 1.0) * sqrt(1.0 - x * x);
  }

  if (l == m) {
    return P_mm;
  }

  double P_prev2 = P_mm;
  double P_prev1 = (2.0 * m + 1.0) * x * P_mm;
  if (l == m + 1) {
    return P_prev1;
  }

  double P_cur = 0.0;
  for (int i = m + 2; i <= l; i++) {
    P_cur = ((2.0 * i - 1.0) * x * P_prev1 - (i + m - 1.0) * P_prev2) / (i - m);
    P_prev2 = P_prev1;
    P_prev1 = P_cur;
  }

  return P_cur;
}

/*
 * Evaluate the plain (unassociated, i.e. m=0) Legendre polynomial P_l at n
 * points x[0..n-1], writing results into P[0..n-1].
 */
void legendre_array(int l, const double *x, int n, double *P) {
  if (!x || !P || n <= 0) {
    return;
  }

  for (int i = 0; i < n; i++) {
    P[i] = legendre(l, x[i]);
  }
}
