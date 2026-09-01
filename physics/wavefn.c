/*
Wave functions, normalization, expectation values
*/
#include "wavefn.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/fft/fft.h" // for momentum expectation
#include "../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double TWO = 2.0;
static const double NORM_TOL = 1e-30;

// Allocation
wavefunction_t *wavefunction_alloc(int n) {
  wavefunction_t *wavefn = malloc(sizeof(wavefunction_t));
  if (!wavefn) {
    return NULL;
  }

  wavefn->psi = cvector_alloc(n);
  if (!wavefn->psi) {
    free(wavefn);

    return NULL;
  }

  wavefn->x = calloc(n, sizeof(double));
  if (!wavefn->x) {
    cvector_free(wavefn->psi);
    free(wavefn);

    return NULL;
  }

  wavefn->n = n;
  wavefn->dx = 1.0;

  return wavefn;
}

void wavefunction_free(wavefunction_t *wavefn) {
  if (!wavefn) {
    return;
  }

  cvector_free(wavefn->psi);
  free(wavefn->x);
  free(wavefn);
}

wavefunction_t *wavefunction_copy(const wavefunction_t *wavefn) {
  if (!wavefn) {
    return NULL;
  }

  wavefunction_t *copy = wavefunction_alloc(wavefn->n);
  if (!copy) {
    return NULL;
  }

  cvector_t *psi_copy = cvector_copy(wavefn->psi);
  if (!psi_copy) {
    wavefunction_free(copy);

    return NULL;
  }

  cvector_free(copy->psi); // free the newly allocated one

  copy->psi = psi_copy;
  memcpy(copy->x, wavefn->x, wavefn->n * sizeof(double));
  copy->dx = wavefn->dx;

  return copy;
}

// Normalization
void wavefunction_normalize(wavefunction_t *wavefn) {
  if (!wavefn) {
    return;
  }

  double norm_sq = 0.0;
  for (int i = 0; i < wavefn->n; i++) {
    norm_sq += c_abs2(wavefn->psi->data[i]);
  }

  norm_sq *= wavefn->dx;
  if (norm_sq < 1e-30) {
    return;
  }

  double norm = sqrt(norm_sq);
  for (int i = 0; i < wavefn->n; i++) {
    wavefn->psi->data[i] = c_scale(wavefn->psi->data[i], 1.0 / norm);
  }
}

// Probability density
double *wavefunction_prob_density(const wavefunction_t *wavefn) {
  if (!wavefn) {
    return NULL;
  }

  double *rho = malloc(wavefn->n * sizeof(double));
  if (!rho) {
    return NULL;
  }

  for (int i = 0; i < wavefn->n; i++) {
    rho[i] = c_abs2(wavefn->psi->data[i]);
  }

  return rho;
}

double wavefunction_prob_in_interval(const wavefunction_t *wavefn, double a,
                                     double b) {
  if (!wavefn) {
    return 0.0;
  }

  double prob = 0.0;
  for (int i = 0; i < wavefn->n; i++) {
    if (wavefn->x[i] >= a && wavefn->x[i] <= b) {
      prob += c_abs2(wavefn->psi->data[i]);
    }
  }

  return prob * wavefn->dx;
}

// Expectation values
double wavefunction_expect_x(const wavefunction_t *wavefn) {
  if (!wavefn) {
    return 0.0;
  }

  double sum = 0.0;
  for (int i = 0; i < wavefn->n; i++) {
    sum += wavefn->x[i] * c_abs2(wavefn->psi->data[i]);
  }

  return sum * wavefn->dx;
}

double wavefunction_expect_x2(const wavefunction_t *wavefn) {
  if (!wavefn) {
    return 0.0;
  }

  double sum = 0.0;
  for (int i = 0; i < wavefn->n; i++) {
    sum += wavefn->x[i] * wavefn->x[i] * c_abs2(wavefn->psi->data[i]);
  }

  return sum * wavefn->dx;
}

// Momentum expectation via FFT
double wavefunction_expect_p(const wavefunction_t *wavefn) {
  if (!wavefn) {
    return 0.0;
  }

  // Transform to momentum space using FFT
  cvector_t *psi_k = position_to_momentum(wavefn->psi, wavefn->dx);
  if (!psi_k) {
    return 0.0;
  }

  // Compute expectation of p = \hbar * k
  double delta_k = TWO * M_PI / (wavefn->n * wavefn->dx);
  double sum = 0.0;
  for (int i = 0; i < wavefn->n; i++) {
    // natural (unshifted) FFT bin order
    double kval = (i < wavefn->n / 2) ? i * delta_k : (i - wavefn->n) * delta_k;

    sum += kval * c_abs2(psi_k->data[i]);
  }

  sum *= delta_k;

  cvector_free(psi_k);

  return sum; // natural units: \hbar=1, so p=k
}

double wavefunction_expect_p2(const wavefunction_t *wavefn) {
  if (!wavefn) {
    return 0.0;
  }

  cvector_t *psi_k = position_to_momentum(wavefn->psi, wavefn->dx);
  if (!psi_k) {
    return 0.0;
  }

  double delta_k = 2.0 * M_PI / (wavefn->n * wavefn->dx);
  double sum = 0.0;
  for (int i = 0; i < wavefn->n; i++) {
    double kval = (i < wavefn->n / 2) ? i * delta_k : (i - wavefn->n) * delta_k;

    sum += kval* kval * c_abs2(psi_k->data[i]);
  }

  sum *= delta_k;

  cvector_free(psi_k);

  return sum;
}

// Uncertainties
double wavefunction_delta_x(const wavefunction_t *wavefn) {
  double x2 = wavefunction_expect_x2(wavefn);
  double x = wavefunction_expect_x(wavefn);

  return sqrt(x2 - x * x);
}

double wavefunction_delta_p(const wavefunction_t *wavefn) {
  double p2 = wavefunction_expect_p2(wavefn);
  double p = wavefunction_expect_p(wavefn);

  return sqrt(p2 - p * p);
}

double wavefunction_uncertainty_product(const wavefunction_t *wavefn) {
  return wavefunction_delta_x(wavefn) * wavefunction_delta_p(wavefn);
}

// I/O
void wavefunction_save(const wavefunction_t *wavefn, const char *filename) {
  if (!wavefn || !filename) {
    return;
  }

  FILE *f = fopen(filename, "w");
  if (!f) {
    fprintf(stderr, "Cannot open %s\n", filename);
    return;
  }

  fprintf(f, "# x  Re(psi)  Im(psi)  |psi|^2\n");
  for (int i = 0; i < wavefn->n; i++) {
    double r = wavefn->psi->data[i].re;
    double im = wavefn->psi->data[i].im;

    fprintf(f, "% .6e  % .6e  % .6e  % .6e\n", wavefn->x[i], r, im,
            r * r + im * im);
  }

  fclose(f);
}

void wavefunction_save_prob(const wavefunction_t *wavefn, const char *filename) {
  if (!wavefn || !filename) {
    return;
  }

  double *rho = wavefunction_prob_density(wavefn);
  if (!rho) {
    return;
  }

  FILE *f = fopen(filename, "w");
  if (!f) {
    free(rho);

    return;
  }

  fprintf(f, "# x  |\\psi|^2\n");
  for (int i = 0; i < wavefn->n; i++) {
    fprintf(f, "% .6e  % .6e\n", wavefn->x[i], rho[i]);
  }

  fclose(f);

  free(rho);
}
