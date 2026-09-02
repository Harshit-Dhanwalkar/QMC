#ifndef QMC_UTIL_H
#define QMC_UTIL_H

#include "matrix.h"
#include "vector.h"

/* Array/Grid utilities */
double *linspace(double start, double end, int num_points);
double *logspace(double start, double end, int num_points);
int *range(int start, int end);

/* Wave function utilities */
void normalize_wavefunction(cvector_t *psi);
double compute_norm_squared(const cvector_t *psi, double grid_dx);

/* Expectation values */
double expectation_value(const cvector_t *psi, const cvector_t *op_psi,
                         double grid_dx);
double expectation_position(const cvector_t *psi, const double *pos_grid,
                            double grid_dx);
double expectation_position_squared(const cvector_t *psi,
                                    const double *pos_grid, double grid_dx);

/* Momentum space (via FFT interface) */
cvector_t *position_to_momentum(const cvector_t *psi_x, double grid_dx);
double expectation_momentum(const cvector_t *psi_k, const double *k_grid,
                            double grid_dk);

/* Statistics */
double mean(const double *data, int n);
double variance(const double *data, int n);
double std_dev(const double *data, int n);

/* Matrix column extraction */
cvector_t *cvector_from_matrix_column(const cmatrix_t *m, int col);

/* File I/O utilities */
void save_wavefunction(const char *filename, const double *pos_grid,
                       const cvector_t *psi, int num_points);
void save_potential(const char *filename, const double *x, const double *V,
                    int n);
void save_eigenvalues(const char *filename, const double *eigenvals, int n);

#endif
