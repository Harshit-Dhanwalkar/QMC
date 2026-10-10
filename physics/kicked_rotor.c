/*
 * Implementation of quantum kicked rotor and classical standard map
 *
 * Quantum propagation (krotor_step):
 *   One kick period is a split-step between two natural bases,
 *     \psi -> \exp(-i (K/\hbar) \cos(\theta))       [ kick, diagonal in angle ]
 *     \psi -> FFT -> \exp(-i \hbar m^2 / 2) -> IFFT [ free rotation,
 *                                                     diagonal in momentum ]
 *
 *   Angle-space kick phase is applied in place; free-rotation phase is applied
 *   between a forward FFT and its inverse. Both phase arrays are precomputed
 *   once per call, so cost per kick is O(N log N) from FFTs and O(N) from phase
 *   multiplications. State is copied into a scratch buffer only once, at top,
 *   so caller's array is untouched until final state is written back after
 *   `steps` kicks
 *
 * Momentum representation (krotor_momentum_distribution):
 *   A normalized forward FFT of angle-space state gives amplitudes c_m in
 *   momentum basis p = \hbar m, with usual FFT index convention: index j holds
 *   m = j for j < N/2 and m = j - N otherwise
 *   Probability array prob[m + N/2] is written in ascending m so caller's
 *   downstream plotting code sees a contiguous "m = -N/2 .. N/2 - 1" axis
 *
 * Classical ensemble (krotor_classical_init / _step / _p2):
 *   Angles uniform in [0, 2\pi) and momenta uniform in [-\hbar/2, \hbar/2)
 *   (width of quantum m = 0 state), advanced by Chirikov standard map
 *      p'      = p + K \sin(\theta)
 *       \theta' = \theta + p' mod 2\pi
 *   RNG is SplitMix64 so ensemble is bit-for-bit reproducible for a given seed
 *   on any platform
 *
 * <m^2> (krotor_momentum2) is computed from momentum distribution returned
 * above so same numbers that feed plot also feed scalar observable; two are
 * guaranteed consistent by construction.
 */
#include "kicked_rotor.h"
#include "../core/fft/fft.h"
#include "complex.h"
#include "vector.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int good_length(const cvector_t *psi) {
  if (!psi || !psi->data || psi->n < 4) {
    return 0;
  }

  int n = psi->n;

  return (n & (n - 1)) == 0;
}

int krotor_init(cvector_t *psi) {
  if (!good_length(psi)) {
    return -1;
  }

  double a = 1.0 / sqrt((double)psi->n);
  cvector_fill(psi, c_new(a, 0.0));

  return 0;
}

int krotor_step(cvector_t *psi, double k, double hbar, int steps) {
  if (!good_length(psi) || !isfinite(k) || !(hbar > 0.0) || !isfinite(hbar) ||
      steps < 1) {
    return -1;
  }

  int n = psi->n;
  cvector_t *work = cvector_copy(psi);

  if (!work) {
    return -2;
  }

  // free-rotation phases \exp(-i \hbar m^2 / 2) and kick phases, precomputed
  complex_t *free_phase = malloc((size_t)n * sizeof(complex_t));
  complex_t *kick_phase = malloc((size_t)n * sizeof(complex_t));

  if (!free_phase || !kick_phase) {
    free(free_phase);
    free(kick_phase);
    cvector_free(work);

    return -2;
  }

  for (int j = 0; j < n; j++) {
    int m = (j < n / 2) ? j : j - n; /* FFT ordering */
    double fp = -0.5 * hbar * (double)m * (double)m;
    double theta = 2.0 * M_PI * j / n;
    double kp = -(k / hbar) * cos(theta);

    free_phase[j] = c_new(cos(fp), sin(fp));
    kick_phase[j] = c_new(cos(kp), sin(kp));
  }

  for (int s = 0; s < steps; s++) {
    for (int j = 0; j < n; j++) {
      work->data[j] = c_mul(work->data[j], kick_phase[j]);
    }

    fft(work);

    for (int j = 0; j < n; j++) {
      work->data[j] = c_mul(work->data[j], free_phase[j]);
    }

    ifft(work);
  }

  for (int j = 0; j < n; j++) {
    psi->data[j] = work->data[j];
  }

  free(free_phase);
  free(kick_phase);
  cvector_free(work);

  return 0;
}

int krotor_momentum_distribution(const cvector_t *psi, double *prob) {
  if (!good_length(psi) || !prob) {
    return -1;
  }

  int n = psi->n;
  cvector_t *work = cvector_copy(psi);
  if (!work) {
    return -2;
  }

  fft_normalized(work);
  for (int j = 0; j < n; j++) {
    // FFT index j holds m = j (j < n/2) or j - n; store at m + n/2
    int m = (j < n / 2) ? j : j - n;

    prob[m + n / 2] = c_abs2(work->data[j]);
  }

  cvector_free(work);

  return 0;
}

double krotor_norm(const cvector_t *psi) {
  if (!good_length(psi)) {
    return NAN;
  }

  double sum = 0.0;
  for (int j = 0; j < psi->n; j++) {
    sum += c_abs2(psi->data[j]);
  }

  return sum;
}

double krotor_momentum2(const cvector_t *psi) {
  if (!good_length(psi)) {
    return NAN;
  }

  int n = psi->n;
  double *prob = malloc((size_t)n * sizeof(double));
  if (!prob) {
    return NAN;
  }

  if (krotor_momentum_distribution(psi, prob) != 0) {
    free(prob);

    return NAN;
  }

  double sum = 0.0;
  double norm = 0.0;

  for (int i = 0; i < n; i++) {
    double m = (double)i - 0.5 * (double)n;

    sum += m * m * prob[i];
    norm += prob[i];
  }

  free(prob);

  return norm > 0.0 ? sum / norm : NAN;
}

static unsigned long long next_u64(unsigned long long *state) {
  unsigned long long z = (*state += 0x9E3779B97F4A7C15ULL);

  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;

  return z ^ (z >> 31);
}

static double unit(unsigned long long *state) {
  return (double)(next_u64(state) >> 11) * (1.0 / 9007199254740992.0);
}

int krotor_classical_init(double *theta, double *p, int n, double hbar,
                          unsigned long long seed) {
  if (!theta || !p || n < 1 || !isfinite(hbar)) {
    return -1;
  }

  unsigned long long state = seed;

  for (int i = 0; i < n; i++) {
    theta[i] = 2.0 * M_PI * unit(&state);
    p[i] = hbar * (unit(&state) - 0.5);
  }

  return 0;
}

int krotor_classical_step(double *theta, double *p, int n, double k,
                          int steps) {
  if (!theta || !p || n < 1 || !isfinite(k) || steps < 1) {
    return -1;
  }

  for (int s = 0; s < steps; s++) {
    for (int i = 0; i < n; i++) {
      p[i] += k * sin(theta[i]);
      theta[i] = fmod(theta[i] + p[i], 2.0 * M_PI);
      if (theta[i] < 0.0) {
        theta[i] += 2.0 * M_PI;
      }
    }
  }

  return 0;
}

double krotor_classical_p2(const double *p, int n) {
  if (!p || n < 1) {
    return NAN;
  }

  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum += p[i] * p[i];
  }

  return sum / n;
}
