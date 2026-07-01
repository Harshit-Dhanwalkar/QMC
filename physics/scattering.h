#ifndef QMC_SCATTERING_H
#define QMC_SCATTERING_H

#include "../core/matrix.h"
#include "../core/vector.h"
#include "potentials.h"

/* Phase shift for a central potential (partial wave) */
double phase_shift(int l, double k, potential_fn V, void *params, double r_min,
                   double r_max, int N);

/* Born approximation for scattering amplitude.
   Returns complex amplitude f(theta).
*/
complex_t born_amplitude(potential_fn V, void *params, double k, double theta,
                         double r_max, int N);

/* Differential cross section from Born amplitude */
double born_cross_section(complex_t f_theta);

#endif
