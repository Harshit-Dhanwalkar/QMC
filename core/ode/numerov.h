#ifndef QMC_NUMEROV_H
#define QMC_NUMEROV_H

#include "../vector.h"

/*
 * Numerov algorithm for 1D time-independent Schr$\"{o}$dinger equation:
 *   -\hbar^2/(2m) d^2\phi/dx^2 + V(x)\phi = E \cdot \phi
 */

typedef struct {
    double *x;          // Position grid
    double *V;          // Potential array V(x)
    int     n;          // Grid points
    double  dx;         // Grid spacing
    double  hbar_sq_2m; // \hbar^2/(2m) in problem units
} numerov_params_t;

typedef struct {
    double    energy;
    cvector_t *psi;
} numerov_solution_t;

/*
 * numerov_shoot: find eigenstate for given level.
 *
 * E_guess is used only to determine the target level index:
 *   level = round(E_guess - V_min - 0.5), clamped to >= 0.
 */
numerov_solution_t *numerov_shoot(numerov_params_t *params, double E_guess,
                                   double E_tol);

/*
 * numerov_integrate: Numerov integration at given energy.
 *
 * Fills psi with wavefunction by forward Numerov integration.
 * Seed: psi[0]=0 (Dirichlet), psi[1]=1e-8.
 * The result is not normalized; call vec_normalize() afterward.
 */
void numerov_integrate(const numerov_params_t *params, double E,
                        cvector_t *psi);

void numerov_solution_free(numerov_solution_t *sol);

#endif
