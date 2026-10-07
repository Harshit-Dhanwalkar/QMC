/*
Quantum walk: coin-and-shift discrete-time quantum walk with a tunable coin
angle, observables, an exact classical binomial reference and asymptotic
variance 1 - \sin(\theta)
*/

#include "quantum_walk.h"
#include "complex.h"
#include "vector.h"
#include <math.h>
#include <stdlib.h>

static int valid(const cvector_t *up, const cvector_t *down) {
  return up && down && up->data && down->data && up->n == down->n && up->n >= 4;
}

int qwalk_init(cvector_t *up, cvector_t *down, int coin) {
  if (!valid(up, down) || coin < QWALK_COIN_UP || coin > QWALK_COIN_SYM) {
    return -1;
  }

  cvector_fill(up, c_new(0.0, 0.0));
  cvector_fill(down, c_new(0.0, 0.0));
  int mid = up->n / 2;

  if (coin == QWALK_COIN_UP) {
    up->data[mid] = c_new(1.0, 0.0);
  } else if (coin == QWALK_COIN_DOWN) {
    down->data[mid] = c_new(1.0, 0.0);
  } else {
    up->data[mid] = c_new(M_SQRT1_2, 0.0);
    down->data[mid] = c_new(0.0, M_SQRT1_2);
  }

  return 0;
}

int qwalk_step(cvector_t *up, cvector_t *down, double theta, int steps) {
  if (!valid(up, down) || !isfinite(theta) || steps < 1) {
    return -1;
  }

  int n = up->n;
  complex_t *tmp = malloc(2 * (size_t)n * sizeof(complex_t));
  if (!tmp) {
    return -2;
  }

  complex_t *ua = tmp;
  complex_t *db = tmp + n;
  double c = cos(theta);
  double s = sin(theta);

  for (int k = 0; k < steps; k++) {
    for (int j = 0; j < n; j++) {
      complex_t a = up->data[j];
      complex_t b = down->data[j];
      // coin, then shift: up moves to j + 1, down moves to j - 1
      int jr = (j + 1) % n;
      int jl = (j + n - 1) % n;

      ua[jr] = c_add(c_scale(a, c), c_scale(b, s));
      db[jl] = c_sub(c_scale(a, s), c_scale(b, c));
    }

    for (int j = 0; j < n; j++) {
      up->data[j] = ua[j];
      down->data[j] = db[j];
    }
  }

  free(tmp);

  return 0;
}

double qwalk_norm(const cvector_t *up, const cvector_t *down) {
  if (!valid(up, down)) {
    return NAN;
  }

  double sum = 0.0;
  for (int j = 0; j < up->n; j++) {
    sum += c_abs2(up->data[j]) + c_abs2(down->data[j]);
  }

  return sum;
}

int qwalk_probability(const cvector_t *up, const cvector_t *down,
                      double *prob) {
  if (!valid(up, down) || !prob) {
    return -1;
  }

  for (int j = 0; j < up->n; j++) {
    prob[j] = c_abs2(up->data[j]) + c_abs2(down->data[j]);
  }

  return 0;
}

double qwalk_mean(const cvector_t *up, const cvector_t *down) {
  double norm = qwalk_norm(up, down);
  if (!(norm > 0.0)) {
    return NAN;
  }

  double sum = 0.0;
  int mid = up->n / 2;

  for (int j = 0; j < up->n; j++) {
    sum += (j - mid) * (c_abs2(up->data[j]) + c_abs2(down->data[j]));
  }

  return sum / norm;
}

double qwalk_variance(const cvector_t *up, const cvector_t *down) {
  double norm = qwalk_norm(up, down);
  if (!(norm > 0.0)) {
    return NAN;
  }

  double m1 = 0.0;
  double m2 = 0.0;
  int mid = up->n / 2;

  for (int j = 0; j < up->n; j++) {
    double p = c_abs2(up->data[j]) + c_abs2(down->data[j]);
    double x = j - mid;

    m1 += x * p;
    m2 += x * x * p;
  }

  m1 /= norm;
  m2 /= norm;

  return m2 - m1 * m1;
}

int qwalk_classical(double *prob, int n, int t) {
  if (!prob || n < 4 || t < 0 || t >= n / 2) {
    return -1;
  }

  for (int j = 0; j < n; j++) {
    prob[j] = 0.0;
  }

  prob[n / 2] = 1.0;
  for (int s = 0; s < t; s++) {
    // in-place Pascal update, sweeping so each source is read before written
    double prev = 0.0;

    for (int j = 0; j < n; j++) {
      double cur = prob[j];
      double next = (j + 1 < n) ? prob[j + 1] : 0.0;

      prob[j] = 0.5 * (prev + next);
      prev = cur;
    }
  }

  return 0;
}

double qwalk_asymptotic_variance(double theta) {
  if (!isfinite(theta)) {
    return NAN;
  }

  return 1.0 - sin(theta);
}
