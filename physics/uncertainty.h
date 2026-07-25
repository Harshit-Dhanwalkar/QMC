#ifndef QMC_UNCERTAINTY_H
#define QMC_UNCERTAINTY_H

#include "potentials.h"
#include "wavefn.h"

/* Container for uncertainty results */
typedef struct {
  double mean_x;
  double mean_p;
  double delta_x;
  double delta_p;
  double product; /* \Delta x * \Delta p */
} uncertainty_t;

/* Compute uncertainties for wavefunction */
uncertainty_t compute_uncertainties(const wavefunction_t *wf);

/* Energy expectation <H> for given potential and particle mass */
double compute_energy_expectation(const wavefunction_t *wf, potential_fn V,
                                  void *params, double mass);

#endif
