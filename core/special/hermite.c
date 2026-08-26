/*
Hermite polynomials (HO eigenstates)
*/
#include "../linalg/tridiag_eigh.h"
#include "../matrix.h"
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

void hermite_array(int n, const double *x, int m, double *H) {
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

/*
 * Zeros of physicists' Hermite polynomial H_n, via Golub-Welsch algorithm
 *
 * NOTE: Ordering: index 0 is largest (most positive) root, decreasing with k
 */
int hermite_zeros_all(int n, double *zeros) {
  if (n <= 0 || !zeros) {
    return 0;
  }
  if (n == 1) {
    zeros[0] = 0.0;
    return 1;
  }

  double *diag = calloc((size_t)n, sizeof(double)); // all zero
  double *offdiag = malloc((size_t)(n - 1) * sizeof(double));
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);

    return 0;
  }

  for (int i = 1; i < n; i++) {
    offdiag[i - 1] = sqrt(i / 2.0);
  }

  eigen_t *eig = tridiag_eigvals(diag, offdiag, n);

  free(diag);
  free(offdiag);
  if (!eig) {
    return 0;
  }

  /* tridiag_eigvals returns ascending order; reverse to match descending
   * (largest-first) convention. */
  for (int k = 0; k < n; k++) {
    zeros[k] = eig->eigenvalues[n - 1 - k];
  }

  eigen_free(eig);

  return n;
}

double hermite_zeros(int n, int k) {
  if (n <= 0 || k < 0 || k >= n) {
    return 0.0;
  }

  double *zeros = malloc((size_t)n * sizeof(double));
  if (!zeros) {
    return 0.0;
  }

  hermite_zeros_all(n, zeros);

  double result = zeros[k];

  free(zeros);

  return result;
}
