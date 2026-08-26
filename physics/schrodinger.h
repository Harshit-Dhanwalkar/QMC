#ifndef QMC_SCHRODINGER_H
#define QMC_SCHRODINGER_H

#include "../core/matrix.h"
#include "../core/ode/numerov.h"
#include "../core/vector.h"
#include "potentials.h"
#include "wavefn.h"

/*
 * Solve 1D TISE by matrix diagonalization (finite difference).
 *
 * Returns eigen_t* with eigenvalues (energies) and eigenvectors
 * (wavefunctions).
 * The eigenvectors are stored as columns of the eigenvectors matrix.
 * The wavefunctions are not normalized to dx; uses wavefunction_normalize.
 */
eigen_t *solve_tise_matrix(double *x, int n, double dx, double hbar_sq_2m,
                           potential_fn V, void *params);

/*
 * Solve TISE by shooting method (Numerov) for a specific bound state.
 *
 * Returns numerov_solution_t with energy and wavefunction.
 * Requires numerov_params_t with grid and potential array.
 */
numerov_solution_t *solve_tise_shoot(numerov_params_t *params, double E_guess,
                                     double E_tol);

/*
 * Solve TISE by Numerov shooting with log-derivative matching at outer
 * classical turning point.
 *
 * NOTE: E_min/E_max must bracket exactly one eigenvalue; n_scan controls how
 * finely that bracket is pre-scanned for a sign change before bisecting.
 *
 * Returns numerov_solution_t with the converged energy and wavefunction, or
 * NULL if no sign change is found in bracket (either no eigenvalue in range, or
 * n_scan too coarse to resolve it).
 */
numerov_solution_t *solve_tise_shoot_matching(numerov_params_t *params,
                                              double E_min, double E_max,
                                              int n_scan, double tol);

/*
 * Time-dependent Schrödinger evolution via Crank-Nicolson.
 *  H is tridiagonal. diag, offdiag are Hamiltonian entries.
 *  \psi is wavefunction (complex) on grid.
 *  dt is time step; steps is number of steps.
 *
 * Returns 0 on success.
 */
int evolve_tdse_crank(const double *diag, const double *offdiag, int n,
                      cvector_t *psi, double dt, int steps);

/*
 * Split-step Fourier method (for potentials that are functions of x only).
 *  Uses FFT to propagate in momentum space.
 *  \psi: input/output wavefunction (complex) on grid.
 *  x: position grid.
 *  V: potential array.
 *  dx: grid spacing.
 *  dt: time step; steps: number of steps.
 *
 * Returns 0 on success.
 */
int evolve_tdse_split_step(cvector_t *psi, const double *x, const double *V,
                           int n, double dx, double dt, int steps, double hbar,
                           double mass);

#endif
