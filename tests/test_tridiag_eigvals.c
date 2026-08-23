/*
 * Test: tridiag_eigvals (core/linalg/tridiag_eigh.c) - eigenvalues-only
 * fast path.
 *
 * 1. Discrete Laplacian (diag=2, offdiag=-1): tridiag_eigvals must match known
 *    closed-form eigenvalues \lambda_k = 2 - 2 * \cos(k * \pi / (n+1)).
 * 2. tridiag_eigvals must agree exactly (same shared core, skipping eigenvector
 *    bookkeeping) with tridiag_eigh's eigenvalues, across several random
 *    tridiagonal matrices.
 * 3. tridiag_eigvals's eigenvectors field must be NULL.
 */

#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.10f expected=%.10f err=%.2e\n", label, got, expected,
         err);

  return err > tol;
}

static int test_against_analytic_laplacian(void) {
  int fail = 0;
  int N = 100;
  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);

  for (int i = 0; i < N; i++) {
    diag[i] = 2.0;
  }

  for (int i = 0; i < N - 1; i++) {
    offdiag[i] = -1.0;
  }

  eigen_t *eig = tridiag_eigvals(diag, offdiag, N);
  if (!eig) {
    printf("  FAIL: tridiag_eigvals returned NULL\n");

    free(diag);
    free(offdiag);

    return 1;
  }

  double max_err = 0.0;
  for (int k = 1; k <= N; k++) {
    double analytic = 2.0 - 2.0 * cos(k * M_PI / (N + 1));
    double err = fabs(eig->eigenvalues[k - 1] - analytic);

    if (err > max_err) {
      max_err = err;
    }
  }

  printf("  max |eigenvalue - analytic| over N=%d: %.3e\n", N, max_err);
  fail |= (max_err > 1e-8);

  fail |= (eig->eigenvectors != NULL);
  if (eig->eigenvectors != NULL) {
    printf("  FAIL: eigenvectors should be NULL for tridiag_eigvals\n");
  }

  eigen_free(eig);
  free(diag);
  free(offdiag);

  return fail;
}

static int test_matches_tridiag_eigh(unsigned int seed, int N) {
  srand(seed);
  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);

  for (int i = 0; i < N; i++) {
    diag[i] = (double)rand() / RAND_MAX * 4.0 - 2.0;
  }

  for (int i = 0; i < N - 1; i++) {
    offdiag[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
  }

  eigen_t *full = tridiag_eigh(diag, offdiag, N);
  eigen_t *vals_only = tridiag_eigvals(diag, offdiag, N);

  int fail = 0;
  if (!full || !vals_only) {
    printf("  FAIL: solver returned NULL (N=%d)\n", N);
    fail = 1;
  } else {
    double max_diff = 0.0;

    for (int i = 0; i < N; i++) {
      double d = fabs(full->eigenvalues[i] - vals_only->eigenvalues[i]);

      if (d > max_diff) {
        max_diff = d;
      }
    }

    char label[64];
    snprintf(label, sizeof label, "N=%d max|tridiag_eigh - tridiag_eigvals|",
             N);
    fail |= check_close(max_diff, 0.0, 1e-10, label);
  }

  if (full) {
    eigen_free(full);
  }

  if (vals_only) {
    eigen_free(vals_only);
  }

  free(diag);
  free(offdiag);

  return fail;
}

int main(void) {
  int failed = 0;

  printf("Discrete Laplacian vs analytic closed form:\n");
  failed += test_against_analytic_laplacian();

  printf("tridiag_eigvals matches tridiag_eigh (random matrices):\n");
  const int Ns[] = {10, 50, 150};
  for (int i = 0; i < 3; i++) {
    failed += test_matches_tridiag_eigh(200 + i, Ns[i]);
  }

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
