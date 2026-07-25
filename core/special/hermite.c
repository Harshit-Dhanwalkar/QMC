/*
Hermite polynomials (HO eigenstates)
*/
#include "special.h"
#include <math.h>
#include <stdlib.h>

// Recurrence: H_0=1, H_1=2x, H_{n+1}=2x H_n - 2n H_{n-1}
double hermite(int n, double x) {
  if (n < 0) {
    return 0.0;
  }
  if (n == 0) {
    return 1.0;
  }
  if (n == 1) {
    return 2.0 * x;
  }

  double H_prev2 = 1.0;
  double H_prev1 = 2.0 * x;
  double H_cur = 0.0;
  for (int i = 1; i < n; i++) {
    H_cur = 2.0 * x * H_prev1 - 2.0 * i * H_prev2;
    H_prev2 = H_prev1;
    H_prev1 = H_cur;
  }

  return H_cur;
}

void hermite_array(int n, double *x, int m, double *H) {
  if (!x || !H || m <= 0) {
    return;
  }

  for (int i = 0; i < m; i++) {
    H[i] = hermite(n, x[i]);
  }
}

// Derivative via recurrence: H'_n = 2n H_{n-1}
double hermite_deriv(int n, double x) {
  if (n <= 0) {
    return 0.0;
  }

  return 2.0 * n * hermite(n - 1, x);
}

// Zeros of H_n (using Newton or asymptotic).
/*  HACK: For simplicity, return approximated zeros.
  For n up to ~100, use asymptotic formula or precomputed values.
  Simple approximation for small n; for larger, use a numerical root finder.
   Implement a simple bisection using recurrence.
   TODO: use Golub-Welsch algorithm
*/
static double hermite_zero_approx(int n, int k) {
  // Asymptotic zeros: for large n, zeros ~\sqrt(2n+1) * \cos(\pi *
  // (4k+3)/(4n+2))
  double N = (double)n;
  double theta = M_PI * (4.0 * k + 3.0) / (4.0 * N + 2.0);

  return sqrt(2.0 * N + 1.0) * cos(theta);
}

double hermite_zeros(int n, int k) {
  if (n <= 0 || k < 0 || k >= n) {
    return 0.0;
  }

  // Use approximation as starting point for Newton
  double x0 = hermite_zero_approx(n, k);
  // Refine with Newton (H_n(x) = 0)
  for (int iter = 0; iter < 10; iter++) {
    double H = hermite(n, x0);
    double Hp = hermite_deriv(n, x0);

    if (fabs(Hp) < 1e-15) {
      break;
    }

    double dx = -H / Hp;
    x0 += dx;
    if (fabs(dx) < 1e-12) {
      break;
    }
  }

  return x0;
}
