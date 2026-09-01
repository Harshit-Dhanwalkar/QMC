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

// static const double GOLDEN_RATIO = (1.0 + sqrt(5.0)) / 2.0;
static const double HALF = 0.5;

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

  return HALF * (a + b);
}

double variational_energy(const wavefunction_t *wavefn, potential_fn pot,
                          void *pot_params, double mass) {
  if (!wavefn || !pot) {
    return 0.0;
  }

  double T = wavefunction_expect_p2(wavefn) / (2.0 * mass);
  double PE = 0.0;

  for (int i = 0; i < wavefn->n; i++) {
    PE += pot(wavefn->x[i], pot_params) * c_abs2(wavefn->psi->data[i]);
  }

  PE *= wavefn->dx;

  return T + PE;
}

static double variational_closure_eval(double alpha, void *closure_ptr) {
  variational_closure_t *c = (variational_closure_t *)closure_ptr;
  c->trial_func(alpha, c->wf);

  return variational_energy(c->wf, c->V, c->params, c->mass);
}

double variational_minimize(double alpha_min, double alpha_max,
                            void (*trial_func)(double alpha,
                                               wavefunction_t *wavefn),
                            wavefunction_t *wavefn, potential_fn pot,
                            void *pot_params, double mass, double tol) {
  if (!trial_func || !wavefn || !pot) {
    return 0.0;
  }

  variational_closure_t closure = {trial_func, wavefn, pot, pot_params, mass};
  double alpha_opt = golden_section_minimize(
      alpha_min, alpha_max, variational_closure_eval, &closure, tol);

  trial_func(alpha_opt, wavefn);

  return variational_energy(wavefn, pot, pot_params, mass);
}
