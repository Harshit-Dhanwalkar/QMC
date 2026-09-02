#ifndef QMC_DRIVEN_H
#define QMC_DRIVEN_H

#include "../core/vector.h"

/*
 * Driven two-level systems: general time-dependent Hamiltonians
 * - Landau-Zener sweeps
 * - Lab-frame driving beyond rotating-wave approximation (RWA).
 *
 * Natural units (\hbar = m = 1)
 */

// Caller-supplied time-dependent scalar (detuning, coupling, etc).
typedef double (*time_fn)(double t, void *params);

// time_fn returning a fixed value
double time_fn_constant(double t, void *params);

// time_fn returning alpha*t (linear sweep)
double time_fn_linear_ramp(double t, void *params);

/*
 * General time-dependent two-level evolution in the rotating-frame form
 * used throughout rabi.c:
 *   i d\psi/dt = H(t) \psi,   H(t) = (1/2) * [[ \Delta(t),  \Omega(t) ],
 *                                             [ \Omega(t), -\Delta(t) ]]
 * integrated via classic RK4. \Delta(t) and \Omega(t) may be any caller-
 * supplied time_fn (detuning sweeps, pulse envelopes, chirps, ...).
 *
 * Where :
 *  \psi: length-2 state vector, evolved in place from t0 to t0 + steps * dt.
 *
 * Returns 0 on success, -1 on invalid input.
 */
typedef struct {
  double t0;
  double dt;
  int steps;
} driven_params_t;

int driven_two_level_evolve(cvector_t *psi, time_fn Delta, void *delta_params,
                            time_fn Omega, void *omega_params,
                            driven_params_t params);

/*
 * Landau-Zener: for a linearly-swept detuning \Delta(t) = alpha * t through a
 * fixed coupling \Omega (an avoided crossing), the closed-form probability of a
 * diabatic transition, i.e. of ending up (as t -> +\infty from t -> -\infty) in
 * bare/diabatic state system started in, rather than adiabatically following
 * instantaneous ground state across crossing is:
 *   P_LZ = \exp( -\pi * \Omega^2 / (2 * \alpha) )
 * NOTE: \alpha > 0 is the sweep rate; \Omega is fixed coupling.
 *   - Large \alpha (fast sweep) -> P_LZ -> 1 (diabatic limit, no time to react)
 *   - Small \alpha (slow sweep) -> P_LZ -> 0 (adiabatic limit, follows
 * instantaneous eigenstate and ends up in other diabatic state).
 */
double landau_zener_probability(double Omega, double alpha);

/*
 * Lab-frame driven two-level system, without rotating-wave approximation:
 *   H(t) = (\omega_0 / 2) * \sigma_z + \Omega_0 * \cos(\omega_L * t + phase) *
 * \sigma_x
 *   - \omega_0: bare transition (angular) frequency
 *   - \Omega_0: physical drive amplitude
 *   - \omega_L: drive (angular) frequency
 *   - phase: drive phase offset.
 * Integrated via RK4.
 *
 * \psi: length-2 state vector, evolved in place from t0 to t0 + steps*dt.
 * Returns 0 on success, -1 on invalid input.
 */
int driven_two_level_evolve_lab_frame(cvector_t *psi, double omega0,
                                      double Omega0, double omega_L,
                                      double phase, double t0, double dt,
                                      int steps);

#endif // QMC_DRIVEN_H
