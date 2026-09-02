#include "utils.h"
#include "complex.h"
#include "fft/fft.h"
#include "matrix.h"
#include "vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

#define LOG_BASE_TEN 10.0
#define PATH_BUFFER_SIZE 512

double *linspace(double start, double end, int num_elements) {
  double *arr = malloc(num_elements * sizeof(double));
  if (!arr) {
    return NULL;
  }

  for (int i = 0; i < num_elements; i++) {
    arr[i] = start + i * (end - start) / (num_elements - 1);
  }

  return arr;
}

double *logspace(double start, double end, int num_elements) {
  double *arr = malloc(num_elements * sizeof(double));
  if (!arr) {
    return NULL;
  }

  for (int i = 0; i < num_elements; i++) {
    double ratio = i / (double)(num_elements - 1);

    arr[i] = pow(LOG_BASE_TEN, start + ratio * (end - start));
  }

  return arr;
}

int *range(int start, int end) {
  int num_elements = end - start;
  int *arr = malloc(num_elements * sizeof(int));
  if (!arr) {
    return NULL;
  }

  for (int i = 0; i < num_elements; i++) {
    arr[i] = start + i;
  }

  return arr;
}

void normalize_wavefunction(cvector_t *psi) {
  double norm = cvector_norm(psi);
  if (norm > 0) {
    for (int i = 0; i < psi->n; i++) {
      psi->data[i] = c_scale(psi->data[i], 1.0 / norm);
    }
  }
}

double compute_norm_squared(const cvector_t *psi, double deltax) {
  double norm_sq = 0.0;
  for (int i = 0; i < psi->n; i++) {
    norm_sq += c_abs2(psi->data[i]);
  }

  return norm_sq * deltax;
}

double expectation_value(const cvector_t *psi, const cvector_t *op_psi,
                         double deltax) {
  if (psi->n != op_psi->n) {
    fprintf(stderr,
            "Error: vector size mismatch in expectation_value\n"); // NOLINT

    return 0.0;
  }

  complex_t result = c_zero();
  for (int i = 0; i < psi->n; i++) {
    result = c_add(result, c_mul(c_conj(psi->data[i]), op_psi->data[i]));
  }

  return result.re * deltax;
}

double expectation_position(const cvector_t *psi, const double *x, double deltax) {
  double expect = 0.0;
  for (int i = 0; i < psi->n; i++) {
    expect += x[i] * c_abs2(psi->data[i]);
  }

  return expect * deltax;
}

double expectation_position_squared(const cvector_t *psi, const double *x,
                                    double deltax) {
  double expect = 0.0;
  for (int i = 0; i < psi->n; i++) {
    expect += x[i] * x[i] * c_abs2(psi->data[i]);
  }

  return expect * deltax;
}

// Momentum expectation: given momentum-space wavefunction \psi_k and k grid
double expectation_momentum(const cvector_t *psi_k, const double *grid_k,
                            double deltak) {
  if (!psi_k || !grid_k) {
    return 0.0;
  }

  double sum = 0.0;
  for (int i = 0; i < psi_k->n; i++) {
    sum += grid_k[i] * c_abs2(psi_k->data[i]);
  }

  return sum * deltak;
}

void save_wavefunction(const char *filename, const double *posx,
                       const cvector_t *psi, int n) {
  char path[PATH_BUFFER_SIZE];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename); // NOLINT

  FILE *file_ptr = fopen(path, "w");
  if (!file_ptr) {
    fprintf(stderr, "Cannot open %s for writing\n", path);

    return;
  }

  fprintf(file_ptr, "# x  |\\psi|^2  Re(\\psi)  Im(\\psi)\n"); // NOLINT
  for (int i = 0; i < n; i++) {
    fprintf(file_ptr, "%.6e  %.6e  %.6e  %.6e\n", posx[i], c_abs2(psi->data[i]),
            psi->data[i].re, psi->data[i].im);
  }

  fclose(file_ptr);
}

void save_eigenvalues(const char *filename, const double *eigenvals,
                      int n_tol) {
  char path[PATH_BUFFER_SIZE];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);

  FILE *file_ptr = fopen(path, "w");
  if (!file_ptr) {
    fprintf(stderr, "Cannot open %s for writing\n", path);

    return;
  }

  fprintf(file_ptr, "# n  Energy\n");
  for (int i = 0; i < n_tol; i++) {
    fprintf(file_ptr, "%d  %.6e\n", i + 1, eigenvals[i]);
  }

  fclose(file_ptr);
}

void save_potential(const char *filename, const double *posx, const double *pot,
                    int n_tol) {
  char path[PATH_BUFFER_SIZE];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);

  FILE *file_ptr = fopen(path, "w");
  if (!file_ptr) {
    fprintf(stderr, "Cannot open %s for writing\n", path);

    return;
  }

  fprintf(file_ptr, "# x  V(x)\n"); // V = pot
  for (int i = 0; i < n_tol; i++) {
    fprintf(file_ptr, "%.6e  %.6e\n", posx[i], pot[i]);
  }

  fclose(file_ptr);
}

/*
 * Transform a position-space wavefunction to momentum space via FFT,
 * normalized to preserve total probability (Parseval's theorem).
 *
 *   \phi(k_j) = (deltax / \sqrt(2 * \pi)) * FFT_raw(\psi_x)[j]
 *
 * using the RAW (unnormalized) forward FFT, not fft_normalized's unitary
 * (1 / \sqrt(N)) convention.
 */

cvector_t *position_to_momentum(const cvector_t *psi_x, double deltax) {
  cvector_t *psi_k = cvector_copy(psi_x);
  if (!psi_k) {
    return NULL;
  }

  fft(psi_k); // raw, unnormalized forward transform

  double scale = deltax / sqrt(2.0 * M_PI);
  for (int i = 0; i < psi_k->n; i++) {
    psi_k->data[i] = c_scale(psi_k->data[i], scale);
  }

  return psi_k;
}

cvector_t *cvector_from_matrix_column(const cmatrix_t *mat, int col) {
  if (!mat || col < 0 || col >= mat->ncols) {
    return NULL;
  }

  cvector_t *vec = cvector_alloc(mat->nrows);
  if (!vec) {
    return NULL;
  }

  for (int i = 0; i < mat->nrows; i++) {
    vec->data[i] = CMAT(mat, i, col);
  }

  return vec;
}
