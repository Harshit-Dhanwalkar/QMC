#ifndef QMC_WAVEFN_H
#define QMC_WAVEFN_H

#include "../core/matrix.h"
#include "../core/vector.h"

typedef struct {
  cvector_t *psi; /* wavefunction values on grid */
  double *x;      /* position grid (owned) */ // NOTE: what is 'owned'?
  double dx;      /* grid spacing */
  int n;          /* number of grid points */
} wavefunction_t;

/* Allocate/free */
wavefunction_t *wavefunction_alloc(int n);
void wavefunction_free(wavefunction_t *wavefn);
wavefunction_t *wavefunction_copy(const wavefunction_t *wavefn);

/* Normalize (\int|\phi|^2 dx = 1) */
void wavefunction_normalize(wavefunction_t *wavefn);

/* Probability density array (|\phi|^2) */
double *wavefunction_prob_density(const wavefunction_t *wavefn);
double wavefunction_prob_in_interval(const wavefunction_t *wavefn, double a,
                                     double b);

/* Expectation values (using central finite differences for p) */
double wavefunction_expect_x(const wavefunction_t *wavefn);
double wavefunction_expect_x2(const wavefunction_t *wavefn);
double wavefunction_expect_p(const wavefunction_t *wavefn); /* uses FFT */
double wavefunction_expect_p2(const wavefunction_t *wavefn);

/* Uncertainty */
double wavefunction_delta_x(const wavefunction_t *wavefn);
double wavefunction_delta_p(const wavefunction_t *wavefn);
double wavefunction_uncertainty_product(const wavefunction_t *wavefn);

/* I/O */
void wavefunction_save(const wavefunction_t *wavefn, const char *filename);
void wavefunction_save_prob(const wavefunction_t *wavefn, const char *filename);

#endif
