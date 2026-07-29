// TODO: implement numerov for potential wall
// numerov_params_t p;
// p.x = x; p.V = V; p.n = n; p.dx = dx;
// p.hbar_sq_2m = HBAR_SQ / (2 * M_ELECTRON);   // in atomic units
//
// numerov_solution_t *sol = numerov_shoot(&p, -0.5, 1e-8);
// if (sol) {
//     printf("Found energy: %f\n", sol->energy);
//     // use sol->psi
//     numerov_solution_free(sol);
// }

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Analytic energy for infinite square well: E_n = n^2 * \pi^2 * \hbar^2/ (2mL^2)
In atomic units (\hbar=m=1): E_n = n^2 \pi^2/ (2L^2)
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
  for (int n = 1; n <= 5; n++)
    printf("    E_%d = %.6f\n", n, E_analytic(n, L));

  // Numerical
  printf("   Numerical Eigenvalue Solution (Finite Difference):\n");

  int M = N - 2; // points
  if (M < 2) {
    free(x);

    return 1;
  }

  double coeff = 1.0 / (2.0 * dx * dx); // \hbar^2 / 2 * m = 0.5; 0.5/dx^2

  cmatrix_t *H = cmatrix_alloc(M, M);
  if (!H) {
    free(x);

    return 1;
  }

  for (int i = 0; i < M; i++) {
    CMAT(H, i, i) = c_real(2.0 * coeff);

    if (i > 0) {
      CMAT(H, i, i - 1) = c_real(-coeff);
    }

    if (i < M - 1) {
      CMAT(H, i, i + 1) = c_real(-coeff);
    }
  }

  eigen_t *eig = cmatrix_eigh(H);
  if (!eig || !eig->eigenvectors) {
    fprintf(stderr, "  FAIL: eigensolver\n");
    cmatrix_free(H);
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
  cmatrix_free(H);
  free(x);

  return 0;
}
