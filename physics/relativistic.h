#ifndef QMC_RELATIVISTIC_H
#define QMC_RELATIVISTIC_H

#include "../core/matrix.h"
#include "../core/vector.h"

/*
 * Solve 1D Klein-Gordon equation for scalar potential V(x), fast/bulk
 * approximate spectrum:
 *   ( -\hbar^2 * c^2 * d^2/dx^2 + m^2 * c^4 ) * \psi = (E - V)^2 * \psi
 */
eigen_t *klein_gordon_1d(double *x, int N, double *V, double m, double hbar,
                         double c);

typedef struct {
  double energy;
  cvector_t *psi; // real-valued (im=0), dx-normalized
  int converged;  // 1 if fixed-point iteration converged within tol
  int iterations;
} klein_gordon_solution_t;

/*
 * Spatially-varying-V(x) Klein-Gordon solve for one target level, via
 * fixed-point iteration:
 *   1. Guess E.
 *   2. Build Linear tridiagonal eigenvalue problem: H(E)\psi = \lambda \psi,
 *       H(E)_diag[i] = 2 * coeff + m^2 * c^4 + 2 * E * V(x_i) - V(x_i)^2
 *   3. Diagonalize, pick eigenvalue closest to previous E^2.
 *   4. E_new = sign(E_guess) * \sqrt(\lambda). Repeat until |E_new-E|<tol or
 *      max_iter is reached.
 *
 * Where
 * x, V    : Grid and potential, as in klein_gordon_1d. N >= 3.
 * E_guess : Initial energy guess; its sign selects energy branch
 *           (positive-energy solutions for E_guess > 0, negative-energy branch
 *           for E_guess < 0) and its magnitude seeds which level gets tracked.
 * tol     : Convergence tolerance on E b/w iterations.
 * max_iter: Iteration cap.
 *
 * Returns NULL on invalid input or allocation failure.
 * NOTE: On success, check sol->converged - iteration may exhaust max_iter
 * without converging for strongly-varying V(x); sol->energy/\psi are
 * last-iterate values regardless.
 */
klein_gordon_solution_t *
klein_gordon_1d_self_consistent(double *x, int N, double *V, double m,
                                double hbar, double c, double E_guess,
                                double tol, int max_iter);

void klein_gordon_solution_free(klein_gordon_solution_t *sol);

/* Solve 1D Dirac equation (2-component) for scalar potential:
   [ c * \sigma_z * p + m * c^2 * \sigma_x + V(x) ] * \phi = E * \phi
   Where \sigma_x, \sigma_z are Pauli matrices.
   Returns eigen_t with eigenvalues (energies) and eigenvectors (spinors).
*/
eigen_t *dirac_1d(double *x, int N, double *V, double m, double hbar, double c);

#endif
