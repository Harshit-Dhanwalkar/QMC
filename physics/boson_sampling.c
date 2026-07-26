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

static double factorial_bs(int n) {
  double r = 1.0;
  for (int i = 2; i <= n; i++) {
    r *= i;
  }

  return r;
}

complex_t boson_sampling_amplitude(const cmatrix_t *U, int M,
                                   const int *input_modes,
                                   const int *output_modes, int N) {
  if (!U || U->nrows != M || U->ncols != M || N < 1 || !input_modes ||
      !output_modes) {
    return c_zero();
  }

  for (int i = 0; i < N; i++) {
    if (input_modes[i] < 0 || input_modes[i] >= M || output_modes[i] < 0 ||
        output_modes[i] >= M) {
      return c_zero();
    }
  }

  // One row-vector orbital per distinct input mode
  cvector_t **mode_rows = calloc(M, sizeof *mode_rows);
  cvector_t **orbitals = malloc(N * sizeof *orbitals);
  if (!mode_rows || !orbitals) {
    free(mode_rows);
    free(orbitals);

    return c_zero();
  }

  int ok = 1;
  for (int i = 0; i < N && ok; i++) {
    int m = input_modes[i];
    if (!mode_rows[m]) {
      mode_rows[m] = cvector_alloc(M);
      if (!mode_rows[m]) {
        ok = 0;

        break;
      }

      for (int j = 0; j < M; j++) {
        mode_rows[m]->data[j] = CMAT(U, m, j);
      }
    }

    orbitals[i] = mode_rows[m];
  }

  complex_t result = c_zero();
  if (ok) {
    complex_t raw = bosonic_permanent_value(orbitals, N, output_modes);

    // Output-mode multiplicity factorials
    double out_mult = 1.0;
    for (int j = 0; j < N; j++) {
      int already_counted = 0;

      for (int k = 0; k < j; k++) {
        if (output_modes[k] == output_modes[j]) {
          already_counted = 1;

          break;
        }
      }

      if (already_counted) {
        continue;
      }

      int count = 1;
      for (int k = j + 1; k < N; k++) {
        if (output_modes[k] == output_modes[j]) {
          count++;
        }
      }

      out_mult *= factorial_bs(count);
    }

    double scale = sqrt(factorial_bs(N) / out_mult);
    result = c_scale(raw, scale);
  }

  for (int m = 0; m < M; m++) {
    cvector_free(mode_rows[m]);
  }

  free(mode_rows);
  free(orbitals);

  return ok ? result : c_zero();
}

double boson_sampling_probability(const cmatrix_t *U, int M,
                                  const int *input_modes,
                                  const int *output_modes, int N) {
  complex_t A = boson_sampling_amplitude(U, M, input_modes, output_modes, N);

  return c_abs2(A);
}

cmatrix_t *beam_splitter_50_50(void) {
  cmatrix_t *U = cmatrix_alloc(2, 2);
  if (!U) {
    return NULL;
  }

  double inv_sqrt2 = 1.0 / sqrt(2.0);
  CMAT(U, 0, 0) = c_real(inv_sqrt2);
  CMAT(U, 0, 1) = c_real(inv_sqrt2);
  CMAT(U, 1, 0) = c_real(inv_sqrt2);
  CMAT(U, 1, 1) = c_real(-inv_sqrt2);

  return U;
}

cmatrix_t *dft_unitary(int M) {
  if (M < 1) {
    return NULL;
  }

  cmatrix_t *U = cmatrix_alloc(M, M);
  if (!U) {
    return NULL;
  }

  double norm = 1.0 / sqrt((double)M);
  for (int j = 0; j < M; j++) {
    for (int k = 0; k < M; k++) {
      double theta = 2.0 * M_PI * j * k / M;
      CMAT(U, j, k) = c_scale(c_new(cos(theta), sin(theta)), norm);
    }
  }

  return U;
}
