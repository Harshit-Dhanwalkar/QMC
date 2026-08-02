/*
 * Non‑degenerate Perturbation Theory
 *
 * Unperturbed: 1D harmonic oscillator (m=1, \hbar=1, \omega=1)
 * Perturbation: V' = \lambda x^4      (anharmonic term)
 * Computes first‑order correction using perturb_nondeg and compares
 * with numerical diagonalisation of full Hamiltonian.
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
  printf(
      " > Non‑degenerate Perturbation Theory (Harmonic + \\lambda x^4)\n\n");

  // Parameters (atomic units: \hbar=m=\omega=1)
  double omega = 1.0;
  double m = 1.0;
  double hbar = 1.0;
  double lambda = 0.05; // perturbation strength

  int N = 501;
  double x_min = -6.0, x_max = 6.0;
  double dx = (x_max - x_min) / (N - 1);
  double *x = linspace(x_min, x_max, N);
  if (!x) {
    fprintf(stderr, "Memory allocation failed\n");

    return 1;
  }

  // Build unperturbed Hamiltonian: H_0 = -(1/2) * d^2/dx^2 + (1/2) * \omega^2 *
  // x^2
  double coeff = -hbar * hbar / (2.0 * m * dx * dx);
  cmatrix_t *H_0 = cmatrix_alloc(N, N);
  if (!H_0) {
    free(x);

    return 1;
  }

  for (int i = 0; i < N; i++) {
    double V0 = 0.5 * m * omega * omega * x[i] * x[i];
    CMAT(H_0, i, i) = c_real(-2.0 * coeff + V0);

    if (i > 0) {
      CMAT(H_0, i, i - 1) = c_real(coeff);
    }

    if (i < N - 1) {
      CMAT(H_0, i, i + 1) = c_real(coeff);
    }
  }

  // Diagonalise H_0 to get unperturbed eigenstates
  eigen_t *eig0 = cmatrix_eigh(H_0);
  if (!eig0) {
    fprintf(stderr, "Eigen decomposition of H_0 failed\n");
    cmatrix_free(H_0);
    free(x);

    return 1;
  }

  // Build perturbation matrix V' = \lambda * x^4 in the basis of H_0 eigenstates 
  // Compute <n|V'|n> by integrating.
  // Compute first-order correction by numerical integration:
  //    E1 = <\phi_n|\lambda x^4|\phi_n>.
  // TODO: Use perturb_nondeg which expects V matrix in eigenbasis.
  // HACK: For simplicity, compute directly.

  // First 5 states and compute <n|x^4|n> using grid
  printf("   State n   E0 (exact)   <x^4>   \\lambda<x^4> (1st order)   E_full "
         "(numerical)   Error (%%)\n");
  printf("   -------   ----------   -------    ----------------   "
         "------------------   --------\n");

  for (int n = 0; n < 5 && n < eig0->n; n++) {
    // Extract eigenvector as cvector
    cvector_t *psi = cvector_from_matrix_column(eig0->eigenvectors, n);
    if (!psi) {
      continue;
    }

    double norm_sq = 0.0;
    for (int i = 0; i < N; i++) {
      norm_sq += c_abs2(psi->data[i]) * dx;
    }

    if (norm_sq > 0.0) {
      double inv = 1.0 / sqrt(norm_sq);
      for (int i = 0; i < N; i++) {
        psi->data[i].re *= inv;
        psi->data[i].im *= inv;
      }
    }

    // Compute <x^4>
    double x4_expect = 0.0;
    for (int i = 0; i < N; i++) {
      double xi = x[i];
      x4_expect += xi * xi * xi * xi * c_abs2(psi->data[i]);
    }
    x4_expect *= dx;

    double E1 = lambda * x4_expect;
    double E0_analytical = hbar * omega * (n + 0.5);

    // Numerical solution: build H = H_0 + \lambda x^4 and diagonalise.
    // For each state by building full H and re-diagonalising.
    // Reuse H_0 and add perturbation.
    cmatrix_t *H_full = cmatrix_copy(H_0);
    if (!H_full) {
      cvector_free(psi);
      continue;
    }

    for (int i = 0; i < N; i++) {
      double Vp = lambda * x[i] * x[i] * x[i] * x[i];
      CMAT(H_full, i, i) = c_add(CMAT(H_full, i, i), c_real(Vp));
    }

    eigen_t *eig_full = cmatrix_eigh(H_full);
    if (!eig_full) {
      cmatrix_free(H_full);
      cvector_free(psi);

      continue;
    }

    double E_full = eig_full->eigenvalues[n];
    double err = fabs(E_full - (E0_analytical + E1)) / fabs(E_full) * 100.0;

    printf(" %2d       %8.6f    %8.6f   %14.6f    %16.6f   %6.2f%%\n", n,
           E0_analytical, x4_expect, E1, E_full, err);

    cmatrix_free(H_full);
    eigen_free(eig_full);
    cvector_free(psi);
  }

  // Cleanup
  cmatrix_free(H_0);
  eigen_free(eig0);
  free(x);

  return 0;
}
