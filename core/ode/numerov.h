#ifndef QMC_NUMEROV_H
#define QMC_NUMEROV_H

#include "../vector.h"

/* Numerov algorithm for 1D time-independent Schrödinger equation:
 * -ℏ**2/(2m) d**2\phi/dx**2 + V(x)\phi = E\cdot\phi
 *
 * Numerov is superior to finite differences for ODE:
 * \phi_{n+1} = (2\phi_n - \phi_{n-1} - (h**2/12)(f_n + 10f_{n-1} +
 * f_{n-2})\phi_n) / (1 - (h**2/12)f_{n+1}) where f_n = (2m/ℏ**2)(E - V_n)
 */

typedef struct {
  double *x;         // Position grid
  double *V;         // Potential array V(x)
  int n;             // Grid points
  double dx;         // Grid spacing
  double hbar_sq_2m; // ℏ**2/(2m) in atomic units
} numerov_params_t;

/* Solve TISE by shooting method
 * Boundary conditions: \phi(x_min) = 0, \phi(x_max) = 0
 * Returns eigenvalue E and eigenfunction \phi
 */
typedef struct {
  double energy;
  cvector_t *psi;
} numerov_solution_t;

numerov_solution_t *numerov_shoot(numerov_params_t *params, double E_guess,
                                  double E_tol);

void numerov_solution_free(numerov_solution_t *sol);

#endif
