#ifndef QMC_RELATIVISTIC_H
#define QMC_RELATIVISTIC_H

#include "../core/matrix.h"
#include "../core/vector.h"
#include "potentials.h"

/*
 * Solve 1D Klein-Gordon equation for scalar potential V(x), fast/bulk
 * approximate spectrum:
 *   ( -\hbar^2 * c^2 * d^2/dx^2 + m^2 * c^4 ) * \psi = (E - V)^2 * \psi
 *
 * NOTE: Approximation: diagonalizes the V-independent (free) operator above,
 * then assigns each level n an energy
 *   V_expect_n + \sqrt(\lambda_n),
 *
 * Where
 * V_expect_n = <n|V|n> is computed from that level's own free-particle
 * eigenvector (1st-order, non-degenerate-perturbation-theory style correction
 * and not single uniform V_avg shift applied to every level).
 *
 * Exact when V(x) is constant (V_expect_n reduces to V_avg identically for
 * every n). For spatially-varying V(x): real improvement over uniform shift for
 * well-separated levels, but being non-degenerate PT can be worse than uniform
 * shift at levels that are quasi-degenerate in free spectrum. Use
 * klein_gordon_1d_self_consistent instead when needed reliable accuracy on V(x)
 * with near-degenerate low-lying levels, or where correctness matters more than
 * speed of single diagonalization.
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

/*
 * Radial Dirac equation for central potential V(r), fixed \kappa
 * (relativistic angular quantum number: \kappa = -(l+1) for j=l+1/2,
 * \kappa = +l for j=l-1/2). Solves coupled first-order radial equations
 * for the "large"/"small" radial components G(r)=r*g(r), F(r)=r*f(r):
 *   \hbar * c * (dG/dr) = -\hbar * c * (\kappa / r ) *G + (E - V + m * c^2) * F
 *   \hbar * c * (dF/dr) =  \hbar * c * (\kappa / r ) *F - (E - V - m * c^2) * G
 * discretized via central differences into a 2N x 2N matrix.
 * Eigenvector columns are 2N-dimensional: rows [0,N) = G(r), rows [N,2N) =
 * F(r).
 *
 * Where
 * r     : uniform radial grid, r[0] > 0 (to avoids kappa/r singularity at
 *         origin)
 * \kappa: relativistic angular quantum number (nonzero integer).
 *
 * Returns eigen_t with 2N eigenvalues/eigenvectors spanning both
 * positive-energy (E > 0) and negative-energy (E < -mc^2-ish) branches, plus
 * grid-truncation artifacts near E=0. Bound states are eigenvalues in (0,
 * m*c^2).
 */
eigen_t *dirac_radial_solve(double *r, int N, int kappa, potential_fn V,
                            void *params, double m, double hbar, double c);

/*
 * Exact relativistic hydrogen/hydrogen-like energy level (Sommerfeld
 * fine-structure formula), for principal quantum number n, relativistic
 * angular quantum number kappa, and nuclear charge Z:
 *   E_{n, \kappa} = m * c^2 * [1 + (Z * \alpha)^2 / (n - |\kappa| +
 * \sqrt(\kappa^2 - (Z * \alpha)^2))^2]^{-1/2}
 *
 * Where
 * \alpha = e_charge^2 / (4 * \pi * eps0* \hbar * c) (i.e. fine-structure
 * constant)
 *
 * Requires n > |\kappa| >= 1 and Z * \alpha < |\kappa| (bound-state regime for
 * point-charge Dirac equation)
 *
 * Returns NAN outside that range.
 */
double dirac_hydrogen_energy_level(int n, int kappa, double Z, double hbar,
                                   double mass, double e_charge, double eps0,
                                   double c);

#endif
