#ifndef QMC_CRANK_NICOLSON_H
#define QMC_CRANK_NICOLSON_H

#include "../complex.h"
#include "../matrix.h"
#include "../vector.h"

/* Crank-Nicolson time evolution for the 1D Schrodinger equation.
   Hamiltonian must be tridiagonal (kinetic + diagonal potential/CAP).
   Solves: (I + i*dt/2*H) \phi_{n+1} = (I - i*dt/2*H) \phi_n
   using a complex tridiagonal Thomas algorithm.
*/

/* Evolve wavefunction psi for one time step dt, time-independent,
   real-valued diagonal H (no CAP).
   diag: diagonal of H (real, size N)
   offdiag: off-diagonal of H (real, size N-1)
   \psi: input/output complex vector
   Returns 0 on success, -1 on failure.
*/
int crank_nicolson_step(const double *diag, const double *offdiag, double dt,
                        cvector_t *psi);

/*
 * General one-step Crank-Nicolson propagation with complex diagonal:
 *   diag[i] = V_diag(x_i) [+ kinetic on-site term] - i * \Gamma(x_i)
 * \Gamma(x_i) >= 0 is an optional complex absorbing potential (CAP); pass
 * an all-zero imaginary part for ordinary (unitary) evolution.
 * offdiag: real kinetic off-diagonal coupling (size N-1) - CAP and any
 * local potential are diagonal-only, so this never needs to be complex.
 * Returns 0 on success, -1 on failure.
 */
int crank_nicolson_step_general(const complex_t *diag, const double *offdiag,
                                double dt, cvector_t *psi);

/* Build tridiagonal Hamiltonian for a 1D TIME-INDEPENDENT system:
   H = -\hbar^2 / (2m) * d^2/dx^2 + V(x)
   discretized: diag[i] = 2*coeff + V[i],
   offdiag[i] = -coeff,
   Where coeff = \hbar^2/(2m) / dx^2.
*/
void build_tridiagonal_hamiltonian(const double *x, const double *V, int N,
                                   double dx, double hbar_sq_2m, double *diag,
                                   double *offdiag);

/* Caller-supplied time-dependent real potential V(x,t). */
typedef double (*potential_time_fn)(double x, double t, void *params);

/*
 * Build complex tridiagonal hamiltonian for time-dependent potential v(x,t),
 * optionally with complex absorbing potential (CAP) layered on top:
 *   diag_out[i]    = 2 * coeff + V(x[i], t, params) - i * absorb[i]
 *   offdiag_out[i] = -coeff
 * absorb: optional (may be NULL) array of length N, \Gamma(x_i) >= 0
 * NULL means no absorption.
 */
void build_tridiagonal_hamiltonian_time_dependent(
    const double *x, int N, double dx, double hbar_sq_2m, potential_time_fn V,
    void *params, double t, const double *absorb, complex_t *diag_out,
    double *offdiag_out);

/*
 * Monomial complex absorbing potential (CAP), Riss-Meyer-style construction:
 * \Gamma ramps smoothly from 0 up to \eta over an absorbing layer of given
 * width at each edge of the grid, and is exactly 0 in interior.
 *   \Gamma(x) = \eta * ((x - x_layer_start)/width)^power  in each layer = 0 in
 * interior width: absolute width (same units as x) of absorbing layer at each
 * edge.
 * \eta: absorption strength at outermost grid point.
 * power: typically 2 or 3
 * absorb_out: length-N array, filled with \Gamma(x_i) >= 0.
 */
void cap_build_monomial(const double *x, int N, double width, double eta,
                        int power, double *absorb_out);

/*
 * Convenience driver: evolve \psi from t0 for `steps` steps of size dt
 * under a time-dependent potential V(x,t) (+ optional CAP absorb, may
 * be NULL), evaluating H at midpoint time of each step for 2nd-order accuracy.
 * \psi is evolved in place.
 * Returns 0 on success, -1 on first failed step.
 */
int crank_nicolson_evolve_time_dependent(const double *x, int N, double dx,
                                         double hbar_sq_2m, potential_time_fn V,
                                         void *params, const double *absorb,
                                         cvector_t *psi, double t0, double dt,
                                         int steps);

#endif
