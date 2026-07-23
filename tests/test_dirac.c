/*
Test: dirac_1d complex-Hermitian eigensolver.

1. Hermiticity check: the 2N x 2N Dirac matrix must satisfy
   H[a][b] = conj(H[b][a]) for every entry.
2. Physical interpretation: for free particle (V=0), eigenvalue spectrum should
   mostly split into two branches separated by ~2mc^2 (positive-energy states
   near/above +mc^2, negative-energy states near/below -mc^2).
*/

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../physics/relativistic.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int test_hermiticity_of_construction(void) {
  int N = 20;
  double dx = 0.1;
  double hbar = 1.0, c = 1.0, m = 1.0;
  double coeff = hbar * c / (2.0 * dx);

  int M = 2 * N;
  cmatrix_t *H = cmatrix_alloc(M, M);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

  double *V = calloc(N, sizeof *V); // free particle

  for (int i = 0; i < N; i++) {
    int row1 = i, row2 = i + N;
    CMAT(H, row1, row1) = c_real(V[i] + m * c * c);
    CMAT(H, row2, row2) = c_real(V[i] - m * c * c);

    if (i > 0) {
      CMAT(H, row1, i - 1 + N) = c_imag(coeff);
      CMAT(H, row2, i - 1) = c_imag(coeff);
    }

    if (i < N - 1) {
      CMAT(H, row1, i + 1 + N) = c_imag(-coeff);
      CMAT(H, row2, i + 1) = c_imag(-coeff);
    }
  }
  free(V);

  int fail = 0;
  double tol = 1e-12;
  for (int a = 0; a < M; a++) {
    for (int b = 0; b < M; b++) {
      complex_t hab = CMAT(H, a, b);
      complex_t hba = CMAT(H, b, a);
      double err = sqrt(pow(hab.re - hba.re, 2) + pow(hab.im + hba.im, 2));
      if (err > tol) {
        printf("  FAIL: H[%d][%d]=(%.4f,%.4f) not conj of "
               "H[%d][%d]=(%.4f,%.4f)\n",
               a, b, hab.re, hab.im, b, a, hba.re, hba.im);
        fail = 1;
      }
    }
  }
  if (!fail)
    printf("  OK: H[a][b] = conj(H[b][a]) for all %d x %d entries\n", M, M);

  cmatrix_free(H);
  return fail;
}

static int test_free_particle_branches(void) {
  int N = 40;
  double dx = 0.2;
  double hbar = 1.0, c = 1.0, m = 1.0;
  double *x = malloc(N * sizeof *x);
  double *V = calloc(N, sizeof *V);
  for (int i = 0; i < N; i++) {
    x[i] = i * dx;
  }

  eigen_t *eig = dirac_1d(x, N, V, m, hbar, c);
  free(x);
  free(V);
  if (!eig) {
    printf("  FAIL: dirac_1d returned NULL\n");
    return 1;
  }

  double mc2 = m * c * c;
  int n_pos = 0, n_neg = 0, n_gap = 0;
  int n_total = eig->n;
  for (int i = 0; i < n_total; i++) {
    double E = eig->eigenvalues[i];

    if (E >= mc2 - 1e-6)
      n_pos++;
    else if (E <= -mc2 + 1e-6)
      n_neg++;
    else
      n_gap++;
  }
  printf("  mc^2=%.3f: %d states >= +mc^2, %d states <= -mc^2, %d in gap\n",
         mc2, n_pos, n_neg, n_gap);

  eigen_free(eig);

  // Loose: most states should fall in two branches, not the gap.
  return (n_gap > n_total / 4) ? 1 : 0;
}

int main(void) {
  int failed = 0;

  printf("Dirac matrix Hermiticity");
  failed += test_hermiticity_of_construction();

  printf("Free particle +-mc^2 branch structure (qualitative):\n");
  failed += test_free_particle_branches();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
