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

double *linspace(double start, double end, int n) {
  double *arr = malloc(n * sizeof(double));
  if (!arr) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    arr[i] = start + i * (end - start) / (n - 1);
  }

  return arr;
}

double *logspace(double start, double end, int n) {
  double *arr = malloc(n * sizeof(double));
  if (!arr) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    double ratio = i / (double)(n - 1);
    arr[i] = pow(10.0, start + ratio * (end - start));
  }

  return arr;
}

int *range(int start, int end) {
  int n = end - start;
  int *arr = malloc(n * sizeof(int));
  if (!arr) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
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

double compute_norm_squared(const cvector_t *psi, double dx) {
  double norm_sq = 0.0;
  for (int i = 0; i < psi->n; i++) {
    norm_sq += c_abs2(psi->data[i]);
  }

  return norm_sq * dx;
}

double expectation_value(const cvector_t *psi, const cvector_t *op_psi,
                         double dx) {
  if (psi->n != op_psi->n) {
    fprintf(stderr, "Error: vector size mismatch in expectation_value\n");

    return 0.0;
  }

  complex_t result = c_zero();
  for (int i = 0; i < psi->n; i++) {
    result = c_add(result, c_mul(c_conj(psi->data[i]), op_psi->data[i]));
  }

  return result.re * dx;
}

double expectation_position(const cvector_t *psi, const double *x, double dx) {
  double expect = 0.0;
  for (int i = 0; i < psi->n; i++) {
    expect += x[i] * c_abs2(psi->data[i]);
  }

  return expect * dx;
}

double expectation_position_squared(const cvector_t *psi, const double *x,
                                    double dx) {
  double expect = 0.0;
  for (int i = 0; i < psi->n; i++) {
    expect += x[i] * x[i] * c_abs2(psi->data[i]);
  }

  return expect * dx;
}

// Momentum expectation: given momentum-space wavefunction \psi_k and k grid
double expectation_momentum(const cvector_t *psi_k, const double *k,
                            double dk) {
  if (!psi_k || !k) {
    return 0.0;
  }

  double sum = 0.0;
  for (int i = 0; i < psi_k->n; i++) {
    sum += k[i] * c_abs2(psi_k->data[i]);
  }

  return sum * dk;
}

void save_wavefunction(const char *filename, const double *x,
                       const cvector_t *psi, int n) {
  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);
  FILE *f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "Cannot open %s for writing\n", path);

    return;
  }

  fprintf(f, "# x  |\\psi|^2  Re(\\psi)  Im(\\psi)\n");
  for (int i = 0; i < n; i++) {
    fprintf(f, "%.6e  %.6e  %.6e  %.6e\n", x[i], c_abs2(psi->data[i]),
            psi->data[i].re, psi->data[i].im);
  }

  fclose(f);
}

void save_eigenvalues(const char *filename, const double *eigenvals, int n) {
  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);
  FILE *f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "Cannot open %s for writing\n", path);

    return;
  }

  fprintf(f, "# n  Energy\n");
  for (int i = 0; i < n; i++) {
    fprintf(f, "%d  %.6e\n", i + 1, eigenvals[i]);
  }

  fclose(f);
}

void save_potential(const char *filename, const double *x, const double *V,
                    int n) {
  char path[512];
  snprintf(path, sizeof path, "%s/%s", QMC_OUTPUT_DIR, filename);
  FILE *f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "Cannot open %s for writing\n", path);

    return;
  }

  fprintf(f, "# x  V(x)\n");
  for (int i = 0; i < n; i++) {
    fprintf(f, "%.6e  %.6e\n", x[i], V[i]);
  }

  fclose(f);
}

cvector_t *position_to_momentum(const cvector_t *psi_x, double dx) {
  cvector_t *psi_k = cvector_copy(psi_x);
  if (!psi_k) {
    return NULL;
  }

  fft_normalized(psi_k);
  fft_shift(psi_k);

  return psi_k;
}

cvector_t *cvector_from_matrix_column(const cmatrix_t *m, int col) {
  if (!m || col < 0 || col >= m->ncols) {
    return NULL;
  }

  cvector_t *v = cvector_alloc(m->nrows);
  if (!v) {
    return NULL;
  }

  for (int i = 0; i < m->nrows; i++) {
    v->data[i] = CMAT(m, i, col);
  }

  return v;
}
