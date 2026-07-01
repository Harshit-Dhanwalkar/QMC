/*
Variational method: compute energy and minimize over parameters.
*/

#include "variational.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "potentials.h"
#include "wavefn.h"
#include <math.h>
#include <stdlib.h>

double variational_energy(const wavefunction_t *wf, potential_fn V,
                          void *params) {
  if (!wf || !V)
    return 0.0;
  // Assume wavefunction is normalized. Compute <H> = <T> + <V>
  double T = wavefunction_expect_p2(wf) / (2.0 * M_ELECTRON);
  double PE = 0.0;
  for (int i = 0; i < wf->n; i++) {
    PE += V(wf->x[i], params) * c_abs2(wf->psi->data[i]);
  }
  PE *= wf->dx;
  return T + PE;
}

double variational_minimize(double alpha_min, double alpha_max,
                            void (*trial_func)(double alpha,
                                               wavefunction_t *wf),
                            wavefunction_t *wf, potential_fn V, void *params,
                            double tol) {
  if (!trial_func || !wf || !V)
    return 0.0;
  // Simple golden section search
  const double phi = (1.0 + sqrt(5.0)) / 2.0;
  double a = alpha_min, b = alpha_max;
  double c = b - (b - a) / phi;
  double d = a + (b - a) / phi;
  trial_func(c, wf);
  double fc = variational_energy(wf, V, params);
  trial_func(d, wf);
  double fd = variational_energy(wf, V, params);
  while (fabs(b - a) > tol) {
    if (fc < fd) {
      b = d;
      d = c;
      c = b - (b - a) / phi;
      trial_func(c, wf);
      fc = variational_energy(wf, V, params);
    } else {
      a = c;
      c = d;
      d = a + (b - a) / phi;
      trial_func(d, wf);
      fd = variational_energy(wf, V, params);
    }
  }
  double alpha_opt = (a + b) / 2.0;
  trial_func(alpha_opt, wf);
  return variational_energy(wf, V, params);
}
