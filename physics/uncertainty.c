/*
\Delta(x)\delta(p), \delta(E)\delta(t), generalized
*/
#include "uncertainty.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "potentials.h"
#include "wavefn.h"

uncertainty_t compute_uncertainties(const wavefunction_t *wf) {
  uncertainty_t u = {0.0, 0.0, 0.0, 0.0, 0.0};
  if (!wf) {
    return u;
  }

  u.mean_x = wavefunction_expect_x(wf);
  u.mean_p = wavefunction_expect_p(wf);
  u.delta_x = wavefunction_delta_x(wf);
  u.delta_p = wavefunction_delta_p(wf);
  u.product = u.delta_x * u.delta_p;
  return u;
}

double compute_energy_expectation(const wavefunction_t *wf, potential_fn V,
                                  void *params, double mass) {
  if (!wf || !V) {
    return 0.0;
  }
  double KE = wavefunction_expect_p2(wf) / (2.0 * mass);
  double PE = 0.0;

  for (int i = 0; i < wf->n; i++) {
    double V_i = V(wf->x[i], params);
    PE += V_i * c_abs2(wf->psi->data[i]);
  }

  PE *= wf->dx;

  return KE + PE;
}
