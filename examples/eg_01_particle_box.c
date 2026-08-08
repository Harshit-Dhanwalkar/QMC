/*
 * Particle in a 1D Infinite Square Well
 */

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Analytic energy for infinite square well:
 *  E_n = n^2 * \pi^2 * \hbar^2 / (2 * m * L^2)
 * In atomic units (\hbar = m =1):
 *  E_n = n^2 * \pi^2 / (2 * L^2)
 */
static double E_analytic(int n, double L) {
  return (double)(n * n) * M_PI * M_PI / (2.0 * L * L);
}

int main(void) {
  printf(" > Particle in a 1D Infinite Square Well\n");

  int N = 51;
  double L = 1.0;
  double dx = L / (N - 1);

  double *x = linspace(0.0, L, N);
  if (!x) {
    return 1;
  }

  // Analytic
  printf("   Analytic energies (atomic units, \\hbar=m=1):\n");
  for (int n = 1; n <= 5; n++) {
    printf("    E_%d = %.6f\n", n, E_analytic(n, L));
  }

  // Numerical: tridiagonal finite-difference
  int M = N - 2; // interior points
  if (M < 2) {
    free(x);

    return 1;
  }

  double coeff = 1.0 / (2.0 * dx * dx); // \hbar^2/(2m) = 0.5 in a.u.
  double *diag = malloc(M * sizeof(double));
  double *offdiag = malloc((M - 1) * sizeof(double));
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);
    free(x);

    return 1;
  }

  for (int i = 0; i < M; i++) {
    diag[i] = 2.0 * coeff;

    if (i < M - 1) {
      offdiag[i] = -coeff;
    }
  }

  eigen_t *eig = tridiag_eigh(diag, offdiag, M);
  free(diag);
  free(offdiag);

  if (!eig) {
    fprintf(stderr, "  FAIL: eigensolver\n");

    free(x);

    return 1;
  }

  printf("  %6s  %12s  %12s  %10s\n", "n", "Numerical", "Analytic", "Error");
  int n_show = (eig->n < 5) ? eig->n : 5;
  for (int i = 0; i < n_show; i++) {
    double E_num = eig->eigenvalues[i];
    double E_ana = E_analytic(i + 1, L);

    printf("  %6d  %12.6f  %12.6f  %10.2e\n", i + 1, E_num, E_ana,
           fabs(E_num - E_ana));
  }

  eigen_free(eig);

  free(x);

  return 0;
}
