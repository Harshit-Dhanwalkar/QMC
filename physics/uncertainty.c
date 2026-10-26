/*
 * \Delta(x)\delta(p), \delta(E)\delta(t), generalized
 */

#include "uncertainty.h"
#include "potentials.h"
#include "variational.h"
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

// Energy expectation value for a given potential V(x) and particle mass m:
//
//   <H> = <T> + <V>
//       = <p^2> / (2m) + \int |\psi(x)|^2 V(x) dx
//
// NOTE: Kinetic term <p^2> is evaluated in momentum space (via FFT) and
// potential term via direct grid quadrature
//
// WARN: This wrapper provides a entry point for energy expectations within
// uncertainty module (e.g., for future \Delta E \Delta t calculations), keeping
// API cohesive for anyone who includes only uncertainty.h
double compute_energy_expectation(const wavefunction_t *wf, potential_fn V,
                                  void *params, double mass) {
  return variational_energy(wf, V, params, mass);
}
