/*
 * Test perturbation theory: compare first-order correction for harmonic + x^4
 * with known analytic result for the ground state.
 * For m= \omega=1, <0|x^4|0> = 3/4.
 */

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../physics/perturbation.h"
#include "../physics/potentials.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Testing perturbation theory (harmonic + \\lambda x^4)\n");

  int N = 501; // grid size
  double x_min = -6.0, x_max = 6.0;
  // TEST:
  // int N = 1501; // grid size
  // double x_min = -12.0, x_max = 12.0;
  double dx = (x_max - x_min) / (N - 1);

  double *x = linspace(x_min, x_max, N);
  if (!x) {
    fprintf(stderr, "FAIL: linspace\n");
    return 1;
  }

  // Build H_0 (finite-difference harmonic oscillator)
  //  H = -\hbar^2/2m \nabla^2 + 1/2 * x^2,
  //  \hbar=m=1 -> coeff = -1 / (2 * dx^2)
  double coeff = -0.5 / (dx * dx);

  cmatrix_t *H_0 = cmatrix_alloc(N, N);
  if (!H_0) {
    fprintf(stderr, "FAIL: cmatrix_alloc\n");
    free(x);

    return 1;
  }

  for (int i = 0; i < N; i++) {
    double V0 = 0.5 * x[i] * x[i];
    // diagonal: kinetic (-2 * coeff > 0) + potential
    CMAT(H_0, i, i) = c_real(-2.0 * coeff + V0);

    if (i > 0) {
      CMAT(H_0, i, i - 1) = c_real(coeff);
    }
    if (i < N - 1) {
      CMAT(H_0, i, i + 1) = c_real(coeff);
    }
  }

  printf("  Diagonalizing %dx%d Hamiltonian...\n", N, N);

  eigen_t *eig = cmatrix_eigh(H_0);
  if (!eig) {
    fprintf(stderr, "FAIL: cmatrix_eigh returned NULL\n");
    cmatrix_free(H_0);
    free(x);

    return 1;
  }

  if (eig->n < 1 || !eig->eigenvectors) {
    fprintf(stderr, "FAIL: eig invalid (n=%d, eigenvectors=%p)\n", eig->n,
            (void *)eig->eigenvectors);
    eigen_free(eig);
    cmatrix_free(H_0);
    free(x);
    return 1;
  }

  // Ground state energy
  printf("    Ground state energy: %.6f (expected ~0.5)\n",
         eig->eigenvalues[0]);
  if (fabs(eig->eigenvalues[0] - 0.5) > 0.05) {
    fprintf(stderr, "FAIL: ground state energy too far from 0.5\n");
    eigen_free(eig);
    cmatrix_free(H_0);
    free(x);

    return 1;
  }

  // Extract ground state \psi_0 (discrete-normalized: \sum|\psi_i|^2 = 1)
  cvector_t *psi0 = cvector_from_matrix_column(eig->eigenvectors, 0);
  if (!psi0) {
    fprintf(stderr, "FAIL: cvector_from_matrix_column returned NULL\n");
    eigen_free(eig);
    cmatrix_free(H_0);
    free(x);

    return 1;
  }

  // Renormalize from discrete (sum = 1) to continuum (sum * dx = 1) convention
  double norm_sq = 0.0;
  for (int i = 0; i < N; i++) {
    norm_sq += c_abs2(psi0->data[i]);
  }
  norm_sq *= dx;

  if (norm_sq < 1e-30) {
    fprintf(stderr, "FAIL: \\psi0 norm is zero\n");
    cvector_free(psi0);
    eigen_free(eig);
    cmatrix_free(H_0);
    free(x);

    return 1;
  }

  double inv_norm = 1.0 / sqrt(norm_sq);
  for (int i = 0; i < N; i++) {
    psi0->data[i].re *= inv_norm;
    psi0->data[i].im *= inv_norm;
  }

  // Compute <x^4> = \int \phi_0^*(x) * x^4 * \phi_0(x) dx
  double x4_expect = 0.0;
  for (int i = 0; i < N; i++) {
    double xi2 = x[i] * x[i];
    x4_expect += xi2 * xi2 * c_abs2(psi0->data[i]);
  }
  x4_expect *= dx;
  printf("    <x^4> = %.6f  (analytic: 0.75)\n", x4_expect);

  // First-order perturbation: \lambda<x^4> comparison
  // For H = H_0 + \lambdax^4, E_0^(1) = \lambda<0|x^4|0> = 0.75 * \lambda
  int passed = (fabs(x4_expect - 0.75) < 0.02);
  printf("    Test %s\n", passed ? "PASSED" : "FAILED");
  if (!passed) {
    printf("    \\Delta = %.6f (tolerance 0.02)\n", fabs(x4_expect - 0.75));
  }

  // Cleanup
  cvector_free(psi0);
  eigen_free(eig);
  cmatrix_free(H_0);
  free(x);

  return passed ? 0 : 1;
}
