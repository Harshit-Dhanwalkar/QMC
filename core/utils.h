#ifndef QMC_UTIL_H
#define QMC_UTIL_H

#include "matrix.h"
#include "vector.h"

/* Array utilities */
double *linspace(double start, double end, int n);
double *logspace(double start, double end, int n);
int *range(int start, int end);

/* Wave function utilities */
void normalize_wavefunction(cvector_t *psi);
double compute_norm_squared(const cvector_t *psi, double dx);

/* Expectation values */
double expectation_value(const cvector_t *psi, const cvector_t *op_psi,
                         double dx);
double expectation_position(const cvector_t *psi, const double *x, double dx);
double expectation_position_squared(const cvector_t *psi, const double *x,
                                    double dx);

/* Momentum space (via FFT interface) */
cvector_t *position_to_momentum(const cvector_t *psi_x, double dx);
double expectation_momentum(const cvector_t *psi_k, const double *k, double dk);

/* Statistics */
double mean(const double *data, int n);
double variance(const double *data, int n);
double std_dev(const double *data, int n);

/* File I/O */
void save_wavefunction(const char *filename, const double *x,
                       const cvector_t *psi, int n);
void save_eigenvalues(const char *filename, const double *eigenvals, int n);
void save_potential(const char *filename, const double *x, const double *V,
                    int n);

cvector_t *cvector_from_matrix_column(const cmatrix_t *m, int col);
#endif
