/*
Laguerre polynomials (Hydrogen)
*/
#include "special.h"
#include <math.h>
#include <stdlib.h>

/* Associated Paguerre P_n^(alpha)(x) using recurrence:
   (n+1)P_{n+1} = (2n + \alpha+1 - x)P_n - (n + \alpha)P_{n-1}
   P_0 = 1, P_1 = \alpha+1 - x.
*/
double laguerre(int n, double alpha, double x) {
  if (n < 0) {
    return 0.0;
  }
  if (n == 0) {
    return 1.0;
  }
  if (n == 1) {
    return (alpha + 1.0) - x;
  }

  double L_prev2 = 1.0;
  double L_prev1 = (alpha + 1.0) - x;
  double L_cur = 0.0;
  for (int i = 1; i < n; i++) {
    L_cur = ((2.0 * i + alpha + 1.0 - x) * L_prev1 - (i + alpha) * L_prev2) /
            (i + 1.0);
    L_prev2 = L_prev1;
    L_prev1 = L_cur;
  }

  return L_cur;
}

void laguerre_array(int n, double alpha, double *x, int N, double *L) {
  if (!x || !L || N <= 0) {
    return;
  }

  for (int i = 0; i < N; i++) {
    L[i] = laguerre(n, alpha, x[i]);
  }
}
