/*
Relativistic QM: Klein-Gordon and Dirac 1D solvers.
*/

#include "relativistic.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/linalg/linalg.h"
#include "../core/matrix.h"
#include <math.h>
#include <stdlib.h>

eigen_t *klein_gordon_1d(double *x, int N, double *V, double m, double hbar,
                         double c) {
  // Discretize: ( -\hbar^2 c^2 d^2/dx^2 + m^2 c^4 ) \psi = (E - V)^2 \psi
  // => H \psi = (E - V)^2 \psi, with H = -\hbar^2 c^2 d^2/dx^2 + m^2 c^4
  double dx = x[1] - x[0];
  double coeff = hbar * hbar * c * c / (dx * dx);
  cmatrix_t *H = cmatrix_alloc(N, N);
  if (!H)
    return NULL;
  for (int i = 0; i < N; i++) {
    CMAT(H, i, i) = c_real(2.0 * coeff + m * m * c * c * c * c);
    if (i > 0)
      CMAT(H, i, i - 1) = c_real(-coeff);
    if (i < N - 1)
      CMAT(H, i, i + 1) = c_real(-coeff);
  }
  // Solve H \psi = \lambda \psi, then E = V + \sqrt(\lambda) (choose positive
  // energy branch)
  eigen_t *eig = cmatrix_eigh_generic(H);
  if (eig) {
    for (int i = 0; i < eig->n; i++) {
      double lambda = eig->eigenvalues[i];
      // Re-evaluate energy: E = V + \sqrt(\lambda) (approximation)
      // NOTE: take average potential as shift.
      double V_avg = 0.0;
      for (int j = 0; j < N; j++)
        V_avg += V[j];
      V_avg /= N;
      eig->eigenvalues[i] = V_avg + sqrt(lambda);
    }
  }
  return eig;
}

eigen_t *dirac_1d(double *x, int N, double *V, double m, double hbar,
                  double c) {
  // Build 2N x 2N matrix:
  // [ V(x) + m c^2,    -i \hbar c d/dx ]
  // [ -i \hbar c d/dx,  V(x) - m c^2   ]
  // using finite differences.
  double dx = x[1] - x[0];
  int M = 2 * N;
  cmatrix_t *H = cmatrix_alloc(M, M);
  if (!H)
    return NULL;
  for (int i = 0; i < M; i++)
    for (int j = 0; j < M; j++)
      CMAT(H, i, j) = c_zero();

  double coeff = -hbar * c / (2.0 * dx);
  for (int i = 0; i < N; i++) {
    int row1 = i, row2 = i + N;
    // Diagonal potential and mass
    CMAT(H, row1, row1) = c_real(V[i] + m * c * c);
    CMAT(H, row2, row2) = c_real(V[i] - m * c * c);
    // Kinetic coupling (off-diagonal)
    if (i > 0) {
      CMAT(H, row1, i - 1 + N) = c_real(coeff);
      CMAT(H, row2, i - 1) = c_real(-coeff);
    }
    if (i < N - 1) {
      CMAT(H, row1, i + 1 + N) = c_real(coeff);
      CMAT(H, row2, i + 1) = c_real(-coeff);
    }
  }
  eigen_t *eig = cmatrix_eigh_generic(H);
  return eig;
}
