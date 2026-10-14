#include "box_revival.h"
#include "complex.h"
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double box_gaussian_coefficients(complex_t *c, int n_modes, double x0,
                                 double sigma, double n0) {
  if (c == NULL || n_modes < 1 || n_modes > BOX_MAX_MODES ||
      !(x0 > 0.0 && x0 < 1.0) || !(sigma >= 0.005 && sigma <= 0.25) ||
      !isfinite(n0)) {
    return -1.0;
  }

  // composite Simpson on [0, 1]; integrand oscillates at most n_modes + n0
  // half-waves
  int steps = 8 * (n_modes + (int)fabs(n0) + 16);
  double h = 1.0 / steps, k0 = M_PI * n0, norm = 0.0;

  steps += steps % 2;
  h = 1.0 / steps;
  for (int n = 1; n <= n_modes; n++) {
    double re = 0.0, im = 0.0;

    for (int i = 0; i <= steps; i++) {
      double x = i * h;
      double wgt = (i == 0 || i == steps) ? 1.0 : (i % 2 ? 4.0 : 2.0);
      double g = exp(-(x - x0) * (x - x0) / (4.0 * sigma * sigma));
      double s = M_SQRT2 * sin(n * M_PI * x) * g;

      re += wgt * s * cos(k0 * x);
      im += wgt * s * sin(k0 * x);
    }

    c[n - 1].re = re * h / 3.0;
    c[n - 1].im = im * h / 3.0;
  }

  // Gaussian is not normalised above: fix it exactly, then report leaked weight
  double gnorm = 0.0;

  for (int i = 0; i <= steps; i++) {
    double x = i * h;
    double wgt = (i == 0 || i == steps) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    double g = exp(-(x - x0) * (x - x0) / (2.0 * sigma * sigma));

    gnorm += wgt * g;
  }

  gnorm *= h / 3.0;
  for (int n = 0; n < n_modes; n++) {
    c[n].re /= sqrt(gnorm);
    c[n].im /= sqrt(gnorm);
    norm += c[n].re * c[n].re + c[n].im * c[n].im;
  }

  double leaked = norm;

  for (int n = 0; n < n_modes; n++) {
    c[n].re /= sqrt(norm);
    c[n].im /= sqrt(norm);
  }

  return leaked;
}

int box_wavefunction(const complex_t *c, int n_modes, double t, const double *x,
                     int nx, complex_t *psi) {
  if (c == NULL || x == NULL || psi == NULL || n_modes < 1 ||
      n_modes > BOX_MAX_MODES || nx < 1 || !isfinite(t)) {
    return -1;
  }

  complex_t ph[BOX_MAX_MODES];

  for (int n = 1; n <= n_modes; n++) {
    // reduce n^2 t modulo 2 \pi through t mod 2 \pi to keep phase accurate for
    // large t
    double tt = fmod(t, BOX_REVIVAL_TIME);

    ph[n - 1] = c_mul(c[n - 1],
                      c_from_polar(1.0, -fmod((double)n * n * tt, 2.0 * M_PI)));
  }

  for (int j = 0; j < nx; j++) {
    double re = 0.0;
    double im = 0.0;

    for (int n = 1; n <= n_modes; n++) {
      double s = M_SQRT2 * sin(n * M_PI * x[j]);

      re += ph[n - 1].re * s;
      im += ph[n - 1].im * s;
    }

    psi[j].re = re;
    psi[j].im = im;
  }

  return 0;
}

int box_density_map(const complex_t *c, int n_modes, const double *x, int nx,
                    double t_max, int nt, double *rho) {
  if (c == NULL || x == NULL || rho == NULL || n_modes < 1 ||
      n_modes > BOX_MAX_MODES || nx < 1 || nt < 2 || !isfinite(t_max)) {
    return -1;
  }

  double *table = malloc(((size_t)n_modes + 2) * (size_t)nx * sizeof *table);

  if (table == NULL) {
    return -1;
  }

  double *re = table + (size_t)n_modes * nx;
  double *im = re + nx;

  for (int n = 1; n <= n_modes; n++) {
    for (int j = 0; j < nx; j++) {
      table[(size_t)(n - 1) * nx + j] = M_SQRT2 * sin(n * M_PI * x[j]);
    }
  }

  for (int i = 0; i < nt; i++) {
    double t = fmod(t_max * i / (nt - 1), BOX_REVIVAL_TIME);

    for (int j = 0; j < nx; j++) {
      re[j] = im[j] = 0.0;
    }

    for (int n = 1; n <= n_modes; n++) {
      if (c[n - 1].re * c[n - 1].re + c[n - 1].im * c[n - 1].im < 1e-16) {
        continue;
      }

      complex_t a = c_mul(
          c[n - 1], c_from_polar(1.0, -fmod((double)n * n * t, 2.0 * M_PI)));
      const double *row = &table[(size_t)(n - 1) * nx];

      for (int j = 0; j < nx; j++) {
        re[j] += a.re * row[j];
        im[j] += a.im * row[j];
      }
    }

    for (int j = 0; j < nx; j++) {
      rho[(size_t)i * nx + j] = re[j] * re[j] + im[j] * im[j];
    }
  }

  free(table);

  return 0;
}

complex_t box_overlap(const complex_t *c, int n_modes, double t) {
  complex_t z = {0.0, 0.0};

  if (c == NULL || n_modes < 1 || n_modes > BOX_MAX_MODES || !isfinite(t)) {
    z.re = z.im = NAN;

    return z;
  }

  double tt = fmod(t, BOX_REVIVAL_TIME);

  for (int n = 1; n <= n_modes; n++) {
    double p = c[n - 1].re * c[n - 1].re + c[n - 1].im * c[n - 1].im,
           a = -fmod((double)n * n * tt, 2.0 * M_PI);

    z.re += p * cos(a);
    z.im += p * sin(a);
  }

  return z;
}

double box_mean_energy(const complex_t *c, int n_modes) {
  if (c == NULL || n_modes < 1 || n_modes > BOX_MAX_MODES) {
    return NAN;
  }

  double e = 0.0;

  for (int n = 1; n <= n_modes; n++) {
    e +=
        (double)n * n * (c[n - 1].re * c[n - 1].re + c[n - 1].im * c[n - 1].im);
  }

  return e;
}

double box_mean_position(const double *x, const complex_t *psi, int nx) {
  if (x == NULL || psi == NULL || nx < 2) {
    return NAN;
  }

  double num = 0.0;
  double den = 0.0;

  for (int j = 0; j + 1 < nx; j++) {
    double dx = x[j + 1] - x[j],
           a = psi[j].re * psi[j].re + psi[j].im * psi[j].im,
           b = psi[j + 1].re * psi[j + 1].re + psi[j + 1].im * psi[j + 1].im;

    num += 0.5 * dx * (a * x[j] + b * x[j + 1]);
    den += 0.5 * dx * (a + b);
  }

  return num / den;
}
