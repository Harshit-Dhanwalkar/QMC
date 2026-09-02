/*
Boson sampling: photon transition amplitudes/probabilities through
linear-optical network
*/

#include "boson_sampling.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "identical.h"
#include <math.h>
#include <stdlib.h>

#define DFT_FACTOR 2.0

static double factorial_bs(int n) {
  double result = 1.0;
  for (int i = 2; i <= n; i++) {
    result *= i;
  }

  return result;
}

static int validate_mode_indices(const int *input_modes,
                                 const int *output_modes, int num_particles,
                                 int num_modes) {
  for (int idx = 0; idx < num_particles; idx++) {
    if (input_modes[idx] < 0 || input_modes[idx] >= num_modes ||
        output_modes[idx] < 0 || output_modes[idx] >= num_modes) {
      return 0;
    }
  }

  return 1;
}

static int build_mode_rows(const cmatrix_t *unitary_matrix,
                           const int *input_modes, int num_particles,
                           cvector_t **mode_rows, cvector_t **orbitals) {
  for (int idx = 0; idx < num_particles; idx++) {
    int mode_idx = input_modes[idx];
    if (!mode_rows[mode_idx]) {
      mode_rows[mode_idx] = cmatrix_get_row(unitary_matrix, mode_idx);
    }

    if (!mode_rows[mode_idx]) {
      return 0;
    }

    orbitals[idx] = mode_rows[mode_idx];
  }

  return 1;
}

static complex_t calculate_permanent(cvector_t **orbitals,
                                     const int *output_modes,
                                     int num_particles) {
  // Logic split to reduce cognitive complexity
  cmatrix_t *submatrix = cmatrix_alloc(num_particles, num_particles);
  if (!submatrix) {
    return c_zero();
  }

  for (int col = 0; col < num_particles; col++) {
    int out_mode = output_modes[col];
    for (int row = 0; row < num_particles; row++) {
      cmatrix_set(submatrix, row, col, cvector_get(orbitals[row], out_mode));
    }
  }

  complex_t perm = cmatrix_permanent(submatrix);
  cmatrix_free(submatrix);

  return perm;
}

complex_t boson_sampling_amplitude(const cmatrix_t *unitary_matrix,
                                   int num_modes, const int *input_modes,
                                   const int *output_modes, int num_particles) {
  if (!unitary_matrix || unitary_matrix->nrows != num_modes ||
      unitary_matrix->ncols != num_modes || num_particles < 1 || !input_modes ||
      !output_modes) {
    return c_zero();
  }

  if (!validate_mode_indices(input_modes, output_modes, num_particles,
                             num_modes)) {
    return c_zero();
  }

  cvector_t **mode_rows =
      (cvector_t **)calloc((size_t)num_modes, sizeof(cvector_t *));
  cvector_t **orbitals =
      (cvector_t **)malloc((size_t)num_particles * sizeof(cvector_t *));

  if (!mode_rows || !orbitals) {
    free((void *)mode_rows);
    free((void *)orbitals);

    return c_zero();
  }

  int is_valid = build_mode_rows(unitary_matrix, input_modes, num_particles,
                                 mode_rows, orbitals);

  complex_t result = c_zero();
  if (is_valid) {
    result = calculate_permanent(orbitals, output_modes, num_particles);

    // Normalize amplitude by occupation factorials
    int *in_counts = (int *)calloc((size_t)num_modes, sizeof(int));
    int *out_counts = (int *)calloc((size_t)num_modes, sizeof(int));
    if (in_counts && out_counts) {
      for (int i = 0; i < num_particles; i++) {
        in_counts[input_modes[i]]++;
        out_counts[output_modes[i]]++;
      }

      double factor = 1.0;
      for (int m = 0; m < num_modes; m++) {
        factor *= factorial_bs(in_counts[m]) * factorial_bs(out_counts[m]);
      }

      result = c_scale(result, 1.0 / sqrt(factor));
    }

    free(in_counts);
    free(out_counts);
  }

  for (int mode_idx = 0; mode_idx < num_modes; mode_idx++) {
    if (mode_rows[mode_idx]) {
      cvector_free(mode_rows[mode_idx]);
    }
  }

  free((void *)mode_rows);
  free((void *)orbitals);

  return is_valid ? result : c_zero();
}

double boson_sampling_probability(const cmatrix_t *unitary_matrix,
                                  int num_modes, const int *input_modes,
                                  const int *output_modes, int num_particles) {
  complex_t amplitude_val = boson_sampling_amplitude(
      unitary_matrix, num_modes, input_modes, output_modes, num_particles);

  return c_abs_sqr(amplitude_val);
}

cmatrix_t *dft_unitary(int num_modes) {
  if (num_modes < 1) {
    return NULL;
  }

  cmatrix_t *unitary_matrix = cmatrix_alloc(num_modes, num_modes);
  if (!unitary_matrix) {
    return NULL;
  }

  double norm = 1.0 / sqrt((double)num_modes);
  for (int row_idx = 0; row_idx < num_modes; row_idx++) {
    for (int col_idx = 0; col_idx < num_modes; col_idx++) {
      double theta = 2.0 * M_PI * row_idx * col_idx / num_modes;

      cmatrix_set(unitary_matrix, row_idx, col_idx,
                  c_scale(c_new(cos(theta), sin(theta)), norm));
    }
  }

  return unitary_matrix;
}

cmatrix_t *beam_splitter_50_50(void) {
  cmatrix_t *unitary_matrix = cmatrix_alloc(2, 2);
  if (!unitary_matrix) {
    return NULL;
  }

  double inv_sqrt2 = 1.0 / sqrt(2.0);
  CMAT(unitary_matrix, 0, 0) = c_real(inv_sqrt2);
  CMAT(unitary_matrix, 0, 1) = c_real(inv_sqrt2);
  CMAT(unitary_matrix, 1, 0) = c_real(inv_sqrt2);
  CMAT(unitary_matrix, 1, 1) = c_real(-inv_sqrt2);

  return unitary_matrix;
}
