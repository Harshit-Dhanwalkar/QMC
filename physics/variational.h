#ifndef QMC_VARIATIONAL_H
#define QMC_VARIATIONAL_H

#include "../core/vector.h"
#include "potentials.h"
#include "wavefn.h"

typedef struct {
  void (*trial_func)(double alpha, wavefunction_t *wf);
  wavefunction_t *wf;
  potential_fn V;
  void *params;
} variational_closure_t;

/*
 * Generic 1D golden-section minimizer of arbitrary scalar function
 * f(x, params) over [a,b].
 */
double golden_section_minimize(double a, double b,
                               double (*f)(double x, void *params),
                               void *params, double tol);

/* Compute expectation value of Hamiltonian for trial wavefunction */
double variational_energy(const wavefunction_t *wf, potential_fn V,
                          void *params);

/*
 * Minimize energy with respect to single parameter alpha
 */
double variational_minimize(double alpha_min, double alpha_max,
                            void (*trial_func)(double alpha,
                                               wavefunction_t *wf),
                            wavefunction_t *wf, potential_fn V, void *params,
                            double tol);

#endif
