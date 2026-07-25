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

// Generic 1D golden-section minimizer
double golden_section_minimize(double a, double b,
                               double (*f)(double x, void *params),
                               void *params, double tol) {
  const double phi = (1.0 + sqrt(5.0)) / 2.0;
  double c = b - (b - a) / phi;
  double d = a + (b - a) / phi;
  double fc = f(c, params);
  double fd = f(d, params);

  while (fabs(b - a) > tol) {
    if (fc < fd) {
      b = d;
      d = c;
      fd = fc;
      c = b - (b - a) / phi;
      fc = f(c, params);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + (b - a) / phi;
      fd = f(d, params);
    }
  }

  return 0.5 * (a + b);
}

double variational_energy(const wavefunction_t *wf, potential_fn V,
                          void *params, double mass) {
  if (!wf || !V) {
    return 0.0;
  }

  double T = wavefunction_expect_p2(wf) / (2.0 * mass);
  double PE = 0.0;

  for (int i = 0; i < wf->n; i++) {
    PE += V(wf->x[i], params) * c_abs2(wf->psi->data[i]);
  }
  PE *= wf->dx;

  return T + PE;
}

static double variational_closure_eval(double alpha, void *closure_ptr) {
  variational_closure_t *c = (variational_closure_t *)closure_ptr;
  c->trial_func(alpha, c->wf);

  return variational_energy(c->wf, c->V, c->params, c->mass);
}

double variational_minimize(double alpha_min, double alpha_max,
                            void (*trial_func)(double alpha,
                                               wavefunction_t *wf),
                            wavefunction_t *wf, potential_fn V, void *params,
                            double mass, double tol) {
  if (!trial_func || !wf || !V) {
    return 0.0;
  }

  variational_closure_t closure = {trial_func, wf, V, params, mass};
  double alpha_opt = golden_section_minimize(
      alpha_min, alpha_max, variational_closure_eval, &closure, tol);

  trial_func(alpha_opt, wf);

  return variational_energy(wf, V, params, mass);
}
