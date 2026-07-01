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

// Allocation
wavefunction_t *wavefunction_alloc(int n) {
  wavefunction_t *wf = malloc(sizeof(wavefunction_t));
  if (!wf)
    return NULL;
  wf->psi = cvector_alloc(n);
  if (!wf->psi) {
    free(wf);
    return NULL;
  }
  wf->x = calloc(n, sizeof(double));
  if (!wf->x) {
    cvector_free(wf->psi);
    free(wf);
    return NULL;
  }
  wf->n = n;
  wf->dx = 1.0;
  return wf;
}

void wavefunction_free(wavefunction_t *wf) {
  if (!wf)
    return;
  cvector_free(wf->psi);
  free(wf->x);
  free(wf);
}

wavefunction_t *wavefunction_copy(const wavefunction_t *wf) {
  if (!wf)
    return NULL;
  wavefunction_t *copy = wavefunction_alloc(wf->n);
  if (!copy)
    return NULL;
  cvector_t *psi_copy = cvector_copy(wf->psi);
  if (!psi_copy) {
    wavefunction_free(copy);
    return NULL;
  }
  cvector_free(copy->psi); // free the newly allocated one
  copy->psi = psi_copy;
  memcpy(copy->x, wf->x, wf->n * sizeof(double));
  copy->dx = wf->dx;
  return copy;
}

// Normalization
void wavefunction_normalize(wavefunction_t *wf) {
  if (!wf)
    return;
  double norm_sq = 0.0;
  for (int i = 0; i < wf->n; i++)
    norm_sq += c_abs2(wf->psi->data[i]);
  norm_sq *= wf->dx;
  if (norm_sq < 1e-30)
    return;
  double norm = sqrt(norm_sq);
  for (int i = 0; i < wf->n; i++)
    wf->psi->data[i] = c_scale(wf->psi->data[i], 1.0 / norm);
}

// Probability density
double *wavefunction_prob_density(const wavefunction_t *wf) {
  if (!wf)
    return NULL;
  double *rho = malloc(wf->n * sizeof(double));
  if (!rho)
    return NULL;
  for (int i = 0; i < wf->n; i++)
    rho[i] = c_abs2(wf->psi->data[i]);
  return rho;
}

double wavefunction_prob_in_interval(const wavefunction_t *wf, double a,
                                     double b) {
  if (!wf)
    return 0.0;
  double prob = 0.0;
  for (int i = 0; i < wf->n; i++) {
    if (wf->x[i] >= a && wf->x[i] <= b)
      prob += c_abs2(wf->psi->data[i]);
  }
  return prob * wf->dx;
}

// Expectation values
double wavefunction_expect_x(const wavefunction_t *wf) {
  if (!wf)
    return 0.0;
  double sum = 0.0;
  for (int i = 0; i < wf->n; i++)
    sum += wf->x[i] * c_abs2(wf->psi->data[i]);
  return sum * wf->dx;
}

double wavefunction_expect_x2(const wavefunction_t *wf) {
  if (!wf)
    return 0.0;
  double sum = 0.0;
  for (int i = 0; i < wf->n; i++)
    sum += wf->x[i] * wf->x[i] * c_abs2(wf->psi->data[i]);
  return sum * wf->dx;
}

// Momentum expectation via FFT
double wavefunction_expect_p(const wavefunction_t *wf) {
  if (!wf)
    return 0.0;
  // Transform to momentum space using FFT
  cvector_t *psi_k = position_to_momentum(wf->psi, wf->dx);
  if (!psi_k)
    return 0.0;
  // Compute expectation of p = \hbar k
  double dk = 2.0 * M_PI / (wf->n * wf->dx);
  double sum = 0.0;
  for (int i = 0; i < wf->n; i++) {
    double k = (i < wf->n / 2) ? i * dk : (i - wf->n) * dk; // after fft_shift
    sum += k * c_abs2(psi_k->data[i]);
  }
  sum *= dk;
  cvector_free(psi_k);
  return sum * HBAR; // HBAR from constants.h
}

double wavefunction_expect_p2(const wavefunction_t *wf) {
  if (!wf)
    return 0.0;
  cvector_t *psi_k = position_to_momentum(wf->psi, wf->dx);
  if (!psi_k)
    return 0.0;
  double dk = 2.0 * M_PI / (wf->n * wf->dx);
  double sum = 0.0;
  for (int i = 0; i < wf->n; i++) {
    double k = (i < wf->n / 2) ? i * dk : (i - wf->n) * dk;
    sum += k * k * c_abs2(psi_k->data[i]);
  }
  sum *= dk;
  cvector_free(psi_k);
  return sum * HBAR * HBAR;
}

// Uncertainties
double wavefunction_delta_x(const wavefunction_t *wf) {
  double x2 = wavefunction_expect_x2(wf);
  double x = wavefunction_expect_x(wf);
  return sqrt(x2 - x * x);
}

double wavefunction_delta_p(const wavefunction_t *wf) {
  double p2 = wavefunction_expect_p2(wf);
  double p = wavefunction_expect_p(wf);
  return sqrt(p2 - p * p);
}

double wavefunction_uncertainty_product(const wavefunction_t *wf) {
  return wavefunction_delta_x(wf) * wavefunction_delta_p(wf);
}

// I/O
void wavefunction_save(const wavefunction_t *wf, const char *filename) {
  if (!wf || !filename)
    return;
  FILE *f = fopen(filename, "w");
  if (!f) {
    fprintf(stderr, "Cannot open %s\n", filename);
    return;
  }
  fprintf(f, "# x  Re(psi)  Im(psi)  |psi|^2\n");
  for (int i = 0; i < wf->n; i++) {
    double r = wf->psi->data[i].re;
    double im = wf->psi->data[i].im;
    fprintf(f, "% .6e  % .6e  % .6e  % .6e\n", wf->x[i], r, im,
            r * r + im * im);
  }
  fclose(f);
}

void wavefunction_save_prob(const wavefunction_t *wf, const char *filename) {
  if (!wf || !filename)
    return;
  double *rho = wavefunction_prob_density(wf);
  if (!rho)
    return;
  FILE *f = fopen(filename, "w");
  if (!f) {
    free(rho);
    return;
  }
  fprintf(f, "# x  |psi|^2\n");
  for (int i = 0; i < wf->n; i++)
    fprintf(f, "% .6e  % .6e\n", wf->x[i], rho[i]);
  fclose(f);
  free(rho);
}
