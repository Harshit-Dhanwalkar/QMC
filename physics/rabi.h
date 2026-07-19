#ifndef QMC_RABI_H
#define QMC_RABI_H

#include "../core/vector.h"

/*
 * Two-level system in the rotating-wave approximation (RWA).
 *
 * H = (\hbar/2) * [[ \Delta, \Omega ], [ \Omega, -\Delta ]]
 * in the {|1> (ground), |2> (excited)} basis. \Delta = \omega_drive -
 * \omega_0 is the detuning (rad/s) of drive from bare transition
 * frequency; \Omega is the real, angular-frequency-units Rabi coupling
 * strength.
 *
 * H is time-independent in rotating frame
 */

// Exact excited-state probability starting from ground state at
// t=0: P_e(t) = (\Omega^2 / \Omega_R^2) * \sin^2(\Omega_R*t/2), where
// \Omega_R = \sqrt(\Omega^2 + \Delta^2) is generalized Rabi frequency.
double rabi_excited_probability(double t, double Omega, double Delta);

/*
 * Exact unitary propagation of two-level spinor \psi (length 2, \psi[0]=ground
 * amplitude, \psi[1]=excited amplitude) from t=0 to time t, using closed-form
 * evolution operator
 *
 * U(t) = \cos(\Omega_R t/2) I - i * \sin(\Omega_R * t/2)/\Omega_R *
 * [[\Delta,\Omega],[\Omega ,-\Delta]]
 *
 * Overwrites psi in place.
 *
 * Returns 0 on success, -1 on invalid input.
 */
int rabi_evolve_exact(cvector_t *psi, double t, double Omega, double Delta);

/*
 * Treating H as a genuine 2x2 real symmetric (trivially tridiagonal) energy
 * Hamiltonian: diag = [\hbar * \Delta/2, -\hbar * \Delta/2], offdiag = [\hbar *
 * \Omega/2]
 */
int rabi_evolve_numerical(cvector_t *psi, double hbar, double Omega,
                          double Delta, double dt, int steps);

#endif
