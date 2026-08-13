/*
 * Test: tridiagonal-solver optimizations in solve_tise_matrix (schrodinger.c)
 * and klein_gordon_1d (relativistic.c).
 *
 * Test proves that by building dense tridiagonal matrix functions
 */

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/linalg/eigen_generic.h"
#include "../core/matrix.h"
#include "../physics/potentials.h"
#include "../physics/relativistic.h"
#include "../physics/schrodinger.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

static int test_solve_tise_matrix_vs_dense(void) {
  int n = 60; // small enough for dense Jacobi to be fast
  double x_min = -5.0, x_max = 5.0;
  double dx = (x_max - x_min) / (n - 1);
  double hbar_sq_2m = 0.5;
  double omega = 1.0;

  double *x = malloc(n * sizeof *x);
  for (int i = 0; i < n; i++) {
    x[i] = x_min + i * dx;
  }

  eigen_t *eig_new =
      solve_tise_matrix(x, n, dx, hbar_sq_2m, V_harmonic, &omega);

  double coeff = hbar_sq_2m / (dx * dx);
  cmatrix_t *H = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    double V_i = V_harmonic(x[i], &omega);
    CMAT(H, i, i) = c_real(2.0 * coeff + V_i);

    if (i > 0) {
      CMAT(H, i, i - 1) = c_real(-coeff);
    }
    if (i < n - 1) {
      CMAT(H, i, i + 1) = c_real(-coeff);
    }
  }
  eigen_t *eig_old = cmatrix_eigh_generic(H);

  cmatrix_free(H);
  free(x);

  if (!eig_new || !eig_old) {
    printf("  FAIL: solver returned NULL\n");
    return 1;
  }

  int fail = 0;
  for (int k = 0; k < 5; k++) {
    char label[32];
    snprintf(label, sizeof label, "eigenvalue[%d]", k);
    fail |= check_close(eig_new->eigenvalues[k], eig_old->eigenvalues[k], 1e-8,
                        label);
  }

  eigen_free(eig_new);
  eigen_free(eig_old);
  return fail;
}

static int test_klein_gordon_vs_dense(void) {
  int N = 60;
  double x_min = -5.0, x_max = 5.0;
  double dx = (x_max - x_min) / (N - 1);
  double m = 1.0, hbar = 1.0, c = 1.0;

  double *x = malloc(N * sizeof *x);
  double *V = calloc(N, sizeof *V); // free particle
  for (int i = 0; i < N; i++) {
    x[i] = x_min + i * dx;
  }

  eigen_t *eig_new = klein_gordon_1d(x, N, V, m, hbar, c);

  double coeff = hbar * hbar * c * c / (dx * dx);
  cmatrix_t *H = cmatrix_alloc(N, N);
  for (int i = 0; i < N; i++) {
    CMAT(H, i, i) = c_real(2.0 * coeff + m * m * c * c * c * c);
    if (i > 0) {
      CMAT(H, i, i - 1) = c_real(-coeff);
    }
    if (i < N - 1) {
      CMAT(H, i, i + 1) = c_real(-coeff);
    }
  }

  eigen_t *eig_old = cmatrix_eigh_generic(H);
  cmatrix_free(H);
  if (eig_old) {
    double V_avg = 0.0;

    for (int j = 0; j < N; j++) {
      V_avg += V[j];
    }

    V_avg /= N;
    for (int i = 0; i < eig_old->n; i++) {
      eig_old->eigenvalues[i] = V_avg + sqrt(eig_old->eigenvalues[i]);
    }
  }

  free(x);
  free(V);

  if (!eig_new || !eig_old) {
    printf("  FAIL: solver returned NULL\n");
    return 1;
  }

  int fail = 0;
  for (int k = 0; k < 5; k++) {
    char label[32];
    snprintf(label, sizeof label, "eigenvalue[%d]", k);
    fail |= check_close(eig_new->eigenvalues[k], eig_old->eigenvalues[k], 1e-6,
                        label);
  }

  eigen_free(eig_new);
  eigen_free(eig_old);

  return fail;
}

int main(void) {
  int failed = 0;

  printf("solve_tise_matrix: optimized vs. dense-Jacobi reference:\n");
  failed += test_solve_tise_matrix_vs_dense();

  printf("klein_gordon_1d: optimized vs. dense-Jacobi reference:\n");
  failed += test_klein_gordon_vs_dense();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
