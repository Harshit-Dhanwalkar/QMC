#include "anderson.h"
#include "complex.h"
#include "vector.h"
#include <math.h>
#include <stdlib.h>

// splitmix64: tiny, deterministic and platform independent
static unsigned long long next_u64(unsigned long long *state) {
  unsigned long long z = (*state += 0x9E3779B97F4A7C15ULL);

  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;

  return z ^ (z >> 31);
}

int anderson_disorder(double *eps, int n, double w, unsigned long long seed) {
  if (!eps || n < 1 || !isfinite(w) || w < 0.0) {
    return -1;
  }

  unsigned long long state = seed;

  for (int j = 0; j < n; j++) {
    // top 53 bits -> uniform in [0, 1)
    double u = (double)(next_u64(&state) >> 11) * (1.0 / 9007199254740992.0);

    eps[j] = w * (u - 0.5);
  }

  return 0;
}

// out = H \psi  (open chain)
static void apply_h(const complex_t *psi, complex_t *out, const double *eps,
                    double hop, int n) {
  for (int j = 0; j < n; j++) {
    complex_t sum = c_scale(psi[j], eps[j]);

    if (j > 0) {
      sum = c_sub(sum, c_scale(psi[j - 1], hop));
    }
    if (j < n - 1) {
      sum = c_sub(sum, c_scale(psi[j + 1], hop));
    }

    out[j] = sum;
  }
}

// out = -i * (H \psi)
static void deriv(const complex_t *psi, complex_t *out, const double *eps,
                  double hop, int n) {
  apply_h(psi, out, eps, hop, n);
  for (int j = 0; j < n; j++) {
    out[j] = c_new(out[j].im, -out[j].re);
  }
}

int anderson_evolve(cvector_t *psi, const double *eps, double hop, double dt,
                    int steps) {
  if (!psi || !psi->data || !eps || psi->n < 2 || !(hop > 0.0) ||
      !isfinite(hop) || !(dt > 0.0) || !isfinite(dt) || steps < 1) {
    return -1;
  }

  int n = psi->n;
  for (int j = 0; j < n; j++) {
    if (!isfinite(eps[j])) {
      return -1;
    }
  }

  size_t len = (size_t)n;
  complex_t *buf = malloc(6 * len * sizeof(complex_t));
  if (!buf) {
    return -2;
  }

  complex_t *k1 = buf;
  complex_t *k2 = buf + len;
  complex_t *k3 = buf + 2 * len;
  complex_t *k4 = buf + 3 * len;
  complex_t *tmp = buf + 4 * len;
  complex_t *cur = buf + 5 * len;

  for (int j = 0; j < n; j++) {
    cur[j] = psi->data[j];
  }

  for (int s = 0; s < steps; s++) {
    deriv(cur, k1, eps, hop, n);
    for (int j = 0; j < n; j++) {
      tmp[j] = c_add(cur[j], c_scale(k1[j], 0.5 * dt));
    }

    deriv(tmp, k2, eps, hop, n);
    for (int j = 0; j < n; j++) {
      tmp[j] = c_add(cur[j], c_scale(k2[j], 0.5 * dt));
    }

    deriv(tmp, k3, eps, hop, n);
    for (int j = 0; j < n; j++) {
      tmp[j] = c_add(cur[j], c_scale(k3[j], dt));
    }

    deriv(tmp, k4, eps, hop, n);
    for (int j = 0; j < n; j++) {
      complex_t incr = c_add(c_add(k1[j], c_scale(k2[j], 2.0)),
                             c_add(c_scale(k3[j], 2.0), k4[j]));

      cur[j] = c_add(cur[j], c_scale(incr, dt / 6.0));
    }
  }

  for (int j = 0; j < n; j++) {
    psi->data[j] = cur[j];
  }

  free(buf);

  return 0;
}

double anderson_norm(const cvector_t *psi) {
  if (!psi || !psi->data) {
    return NAN;
  }

  double sum = 0.0;
  for (int j = 0; j < psi->n; j++) {
    sum += c_abs2(psi->data[j]);
  }

  return sum;
}

double anderson_width(const cvector_t *psi) {
  double norm = anderson_norm(psi);
  if (!(norm > 0.0)) {
    return NAN;
  }

  double m1 = 0.0;
  double m2 = 0.0;

  for (int j = 0; j < psi->n; j++) {
    double p = c_abs2(psi->data[j]);

    m1 += p * j;
    m2 += p * j * j;
  }

  m1 /= norm;
  m2 /= norm;

  double var = m2 - m1 * m1;

  return var > 0.0 ? sqrt(var) : 0.0;
}

double anderson_ipr(const cvector_t *psi) {
  double norm = anderson_norm(psi);
  if (!(norm > 0.0)) {
    return NAN;
  }

  double sum4 = 0.0;

  for (int j = 0; j < psi->n; j++) {
    double p = c_abs2(psi->data[j]);

    sum4 += p * p;
  }

  return sum4 / (norm * norm);
}

double anderson_energy(const cvector_t *psi, const double *eps, double hop) {
  double norm = anderson_norm(psi);
  if (!eps || !(hop > 0.0) || !(norm > 0.0)) {
    return NAN;
  }

  double sum = 0.0;
  int n = psi->n;
  for (int j = 0; j < n; j++) {
    complex_t a = psi->data[j];

    sum += eps[j] * c_abs2(a);
    if (j < n - 1) {
      complex_t b = psi->data[j + 1];

      sum -= 2.0 * hop * (a.re * b.re + a.im * b.im);
    }
  }

  return sum / norm;
}

double anderson_lyapunov(const double *eps, int n, double hop, double energy) {
  if (!eps || n < 2 || !(hop > 0.0) || !isfinite(energy)) {
    return NAN;
  }

  double prev = 0.0;
  double cur = 1.0;
  double acc = 0.0;

  for (int j = 0; j < n; j++) {
    double next = (energy - eps[j]) / hop * cur - prev;

    prev = cur;
    cur = next;

    double scale = sqrt(prev * prev + cur * cur);

    acc += log(scale);
    prev /= scale;
    cur /= scale;
  }

  return acc / n;
}

double anderson_lyapunov_weak(double w, double hop, double energy) {
  double gap = 4.0 * hop * hop - energy * energy;
  if (!isfinite(w) || w < 0.0 || !(hop > 0.0) || !(gap > 0.0)) {
    return NAN;
  }

  return w * w / (24.0 * gap);
}
