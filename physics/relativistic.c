/*
Relativistic QM: Klein-Gordon and Dirac 1D solvers.
*/

#include "relativistic.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/linalg/linalg.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include <math.h>
#include <stdlib.h>

eigen_t *klein_gordon_1d(double *x, int N, double *V, double m, double hbar,
                         double c) {
  if (!x || !V || N < 2) {
    return NULL;
  }

  // Discretize: ( -\hbar^2 c^2 d^2/dx^2 + m^2 c^4 ) \psi = (E - V)^2 \psi
  double dx = x[1] - x[0];
  double coeff = hbar * hbar * c * c / (dx * dx);
  double mc2 = m * c * c;

  double V_avg = 0.0;
  for (int j = 0; j < N; j++) {
    V_avg += V[j];
  }
  V_avg /= N;

  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);

    return NULL;
  }

  for (int i = 0; i < N; i++) {
    diag[i] = 2.0 * coeff + mc2 * mc2;
    if (i < N - 1) {
      offdiag[i] = -coeff;
    }
  }

  eigen_t *eig = tridiag_eigh(diag, offdiag, N);
  free(diag);
  free(offdiag);
  if (!eig) {
    return NULL;
  }

  // E = V_avg + sqrt(lambda) (positive-energy branch)
  for (int i = 0; i < eig->n; i++) {
    eig->eigenvalues[i] = V_avg + sqrt(eig->eigenvalues[i]);
  }

  return eig;
}

klein_gordon_solution_t *
klein_gordon_1d_self_consistent(double *x, int N, double *V, double m,
                                double hbar, double c, double E_guess,
                                double tol, int max_iter) {
  if (!x || !V || N < 3 || max_iter < 1) {
    return NULL;
  }

  double dx = x[1] - x[0];
  double coeff = hbar * hbar * c * c / (dx * dx);
  double mc4 = m * m * c * c * c * c;
  int sign = (E_guess >= 0.0) ? 1 : -1;

  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);

    return NULL;
  }
  for (int i = 0; i < N - 1; i++) {
    offdiag[i] = -coeff;
  }

  double E = E_guess;
  double lambda_target = E * E;
  eigen_t *eig = NULL;
  int level = -1;
  int converged = 0;
  int iter;

  for (iter = 0; iter < max_iter; iter++) {
    // H(E)_diag[i] = 2*coeff + m^2c^4 + V(x_i)^2 - 2*E*V(x_i)
    for (int i = 0; i < N; i++) {
      diag[i] = 2.0 * coeff + mc4 + V[i] * V[i] - 2.0 * E * V[i];
    }

    if (eig) {
      eigen_free(eig);
    }
    eig = tridiag_eigh(diag, offdiag, N);
    if (!eig) {
      free(diag);
      free(offdiag);

      return NULL;
    }

    level = 0;
    double best = fabs(eig->eigenvalues[0] - lambda_target);
    for (int k = 1; k < eig->n; k++) {
      double diff = fabs(eig->eigenvalues[k] - lambda_target);
      if (diff < best) {
        best = diff;
        level = k;
      }
    }

    double lambda = eig->eigenvalues[level];
    if (lambda < 0.0) {
      // No real-energy solution at this iterate for tracked level
      break;
    }

    double E_new = sign * sqrt(lambda);
    if (fabs(E_new - E) < tol) {
      E = E_new;
      converged = 1;
      iter++;

      break;
    }

    E = E_new;
    lambda_target = E * E;
  }

  free(diag);
  free(offdiag);

  if (!eig || level < 0) {
    if (eig) {
      eigen_free(eig);
    }

    return NULL;
  }

  klein_gordon_solution_t *sol = malloc(sizeof *sol);
  if (!sol) {
    eigen_free(eig);

    return NULL;
  }

  cvector_t *psi = cvector_alloc(N);
  if (!psi) {
    free(sol);
    eigen_free(eig);

    return NULL;
  }

  for (int i = 0; i < N; i++) {
    psi->data[i] = CMAT(eig->eigenvectors, i, level);
  }
  eigen_free(eig);

  // dx-weighted normalization
  double norm = 0.0;
  for (int i = 0; i < N; i++) {
    norm += psi->data[i].re * psi->data[i].re;
  }
  norm *= dx;
  if (norm > 1e-300 && isfinite(norm)) {
    double inv = 1.0 / sqrt(norm);
    for (int i = 0; i < N; i++) {
      psi->data[i].re *= inv;
    }
  }

  sol->energy = E;
  sol->psi = psi;
  sol->converged = converged;
  sol->iterations = iter;

  return sol;
}

void klein_gordon_solution_free(klein_gordon_solution_t *sol) {
  if (!sol) {
    return;
  }
  cvector_free(sol->psi);
  free(sol);
}

eigen_t *dirac_1d(double *x, int N, double *V, double m, double hbar,
                  double c) {
  // Build 2N x 2N Hermitian matrix:
  // [ V(x) + mc^2,          -i * \hbar * c d/dx ]
  // [ -i * hbar * c * d/dx,  V(x) - mc^2        ]
  // via central differences.
  double dx = x[1] - x[0];
  int M = 2 * N;
  cmatrix_t *H = cmatrix_alloc(M, M);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

  double coeff = hbar * c / (2.0 * dx);
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

  eigen_t *eig = cmatrix_eigh_complex(H);
  cmatrix_free(H);

  return eig;
}
