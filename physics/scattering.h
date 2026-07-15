#ifndef QMC_SCATTERING_H
#define QMC_SCATTERING_H

#include "../core/complex.h"
#include "potentials.h"

/*
 * phase_shift: scattering phase shift for partial wave l, via Numerov
 * integration of radial equation + asymptotic
 *
 * l:        partial wave (angular momentum)
 * k:        wavenumber, E = hbar_sq_2m * k^2
 * V, params: central potential V(r) and its parameters
 * r_min:    small nonzero starting radius
 * r_max:    extend well beyond range of V
 * N:        number of radial grid points
 * hbar_sq_2m: hbar^2/(2m)
 *
 * Returns the phase shift delta_l via atan2, in (-pi, pi].
 */
double phase_shift(int l, double k, potential_fn V, void *params, double r_min,
                   double r_max, int N, double hbar_sq_2m);

/* Born approximation scattering amplitude f(\theta):
   f(θ) = - (2m / \hbar^2) (1/q) \int_0^\infty r * V(r) * \sin(q*r) dr,
   q = 2k * \sin(\theta / 2)
*/
complex_t born_amplitude(potential_fn V, void *params, double k, double theta,
                         double r_max, int N);

// Differential cross section: d\sigma/d\Omega = |f(\theta)|^2
double born_cross_section(complex_t f_theta);

#endif
