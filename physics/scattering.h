#ifndef QMC_SCATTERING_H
#define QMC_SCATTERING_H

#include "../core/complex.h"
#include "potentials.h"

/*
 * phase_shift: scattering phase shift for partial wave l, via Numerov
 * integration of radial equation + asymptotic matching.
 *
 * l         : partial wave (angular momentum)
 * k         : wavenumber, E = hbar_sq_2m * k^2
 * V, params : central potential V(r) and its parameters
 * r_min     : small nonzero starting radius
 * r_max     : extend well beyond range of V
 * N         : number of radial grid points
 * hbar_sq_2m: \hbar^2/(2m)
 *
 * Returns phase shift \delta_l via \atan2, in (-\pi, \pi].
 */
double phase_shift(int l, double k, potential_fn V, void *params, double r_min,
                   double r_max, int N, double hbar_sq_2m);

/* Born approximation scattering amplitude f(\theta):
   f(\theta) = - (1/hbar_sq_2m) * (1/q) * \int_0^\infty r * V(r) * \sin(q * r)
   dr, q = 2k * \sin(\theta / 2) */
complex_t born_amplitude(potential_fn V, void *params, double k, double theta,
                         double r_max, int N, double hbar_sq_2m);

// Differential cross section: d\sigma / d\Omega = |f(\theta)|^2
double born_cross_section(complex_t f_theta);

/* Partial-wave observables built from an array of phase shifts */

/* S-matrix element for partial wave l: S_l = \exp(2 * i * \delta_l) */
complex_t s_matrix_element(double delta_l);

/* Full scattering amplitude from a partial-wave phase-shift set:
 *   f(\theta) = (1/k) * \sum_{l=0}^{l_max} (2l+1) * \exp^{i * \delta_l} *
 * \sin(\delta_l) P_l(\cos(\theta)) * \delta_l
 * must be an array of length l_max+1 */
complex_t partial_wave_amplitude(const double *delta_l, int l_max, double k,
                                 double theta);

/*
 * Total cross section from a partial-wave phase-shift set:
 *   \sigma_{tot} = (4 * \pi / k^2) * \sum_l (2l+1) \sin^2(\delta_l)
 * Cross-check two against optical theorem: \sigma_{tot} = (4 * \pi / k) *
 * Im[f(0)]
 */
double total_cross_section_partial_waves(const double *delta_l, int l_max,
                                         double k);

/* Lippmann-Schwinger / Born-series scattering
 * These solve radial LS equation for u_l(r) on a dense grid */

/* Exact (non-perturbative) phase shift via direct linear solve of discretized
 * Lippmann-Schwinger equation */
double phase_shift_lippmann_schwinger(int l, double k, potential_fn V,
                                      void *params, double r_min, double r_max,
                                      int N, double hbar_sq_2m);

/*
 * Phase shift via truncated Neumann/Born series to a given order:
 *   order = 1 reduces exactly to standard first-Born partial-wave phase shift
 *           (no iteration, u = hat_j_l(kr)).
 *   order >= 2 self-consistently iterates u through the LS kernel; with more
 *           decimal places.
 *   As order -> \infty this converges to phase_shift_lippmann_schwinger()'s
 *           exact result (for potentials where Born  series converges at all)
 */
double born_series_phase_shift(int l, double k, potential_fn V, void *params,
                               double r_min, double r_max, int N,
                               double hbar_sq_2m, int order);

#endif
