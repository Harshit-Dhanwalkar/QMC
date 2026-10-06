/*
Monte Carlo sampling of hydrogen orbitals
*/

#include "orbital_sample.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/special/special.h"
#include "../core/vector.h"
#include "hydrogen.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { RADIAL_TABLE = 4096, MAX_REJECT = 20000, BOUND_GRID = 96 };

static int valid_quantum_numbers(int n, int l, int m) {
  return n >= 1 && n <= ORBITAL_MAX_N && l >= 0 && l < n && abs(m) <= l;
}

// xorshift64*: small, deterministic, identical everywhere
static double next_uniform(unsigned long long *s) {
  *s ^= *s >> 12;
  *s ^= *s << 25;
  *s ^= *s >> 27;

  return (double)((*s * 2685821657736338717ULL) >> 11) * (1.0 / 9007199254740992.0);
}

// Value of the (complex or real) harmonic.
static complex_t harmonic(int l, int m, int real_form, double theta,
                          double phi) {
  if (!real_form) {
    return spherical_harmonic(l, m, theta, phi);
  }
  if (m == 0) {
    return c_real(spherical_harmonic_real(l, 0, theta, phi));
  }
  if (m > 0) {
    return c_real(sqrt(2.0) * spherical_harmonic_real(l, m, theta, phi));
  }

  return c_real(sqrt(2.0) * spherical_harmonic_imag(l, -m, theta, phi));
}

double hydrogen_orbital_angular_density(int l, int m, int real_form,
                                        double theta, double phi) {
  if (l < 0 || abs(m) > l) {
    return NAN;
  }

  return c_abs2(harmonic(l, m, real_form, theta, phi));
}

double hydrogen_orbital_radial_au(int n, int l, double r) {
  if (n < 1 || l < 0 || l >= n || !(r >= 0.0)) {
    return NAN;
  }

  double r_si = r * AU_LENGTH;
  cvector_t *psi = hydrogen_radial_wavefunction(&r_si, 1, n, l);
  if (!psi) {
    return NAN;
  }
  double value = psi->data[0].re * pow(AU_LENGTH, 1.5);

  cvector_free(psi);

  return value;
}

int hydrogen_orbital_sample(int n, int l, int m, int real_form, int count,
                            unsigned long long seed, double *out) {
  if (!valid_quantum_numbers(n, l, m) || count < 1 || !out) {
    return -1;
  }

  // Radial grid out to where the n-th shell has long decayed
  double rmax = (double)n * (2.0 * n + 16.0);
  double dr = rmax / RADIAL_TABLE;
  double *r_si = malloc(RADIAL_TABLE * sizeof *r_si);
  double *cdf = malloc(RADIAL_TABLE * sizeof *cdf);
  signed char *sign = malloc(RADIAL_TABLE);
  if (!r_si || !cdf || !sign) {
    free(r_si);
    free(cdf);
    free(sign);

    return -2;
  }

  for (int i = 0; i < RADIAL_TABLE; i++) {
    r_si[i] = (i + 0.5) * dr * AU_LENGTH;
  }

  cvector_t *rad = hydrogen_radial_wavefunction(r_si, RADIAL_TABLE, n, l);
  if (!rad) {
    free(r_si);
    free(cdf);
    free(sign);

    return -2;
  }

  double scale = pow(AU_LENGTH, 1.5);
  double total = 0.0;
  for (int i = 0; i < RADIAL_TABLE; i++) {
    double R = rad->data[i].re * scale;
    double r = (i + 0.5) * dr;

    total += R * R * r * r;
    cdf[i] = total;
    sign[i] = (R < 0.0) ? -1 : 1;
  }

  cvector_free(rad);
  free(r_si);

  // Upper bound of |Y|^2 for rejection (grid max with margin)
  double bound = 0.0;
  for (int a = 0; a <= BOUND_GRID; a++) {
    for (int b = 0; b < BOUND_GRID; b++) {
      double th = M_PI * a / BOUND_GRID;
      double ph = 2.0 * M_PI * b / BOUND_GRID;
      double d = c_abs2(harmonic(l, m, real_form, th, ph));

      if (d > bound) {
        bound = d;
      }
    }
  }

  bound *= 1.15;

  unsigned long long state = seed ? seed : 0x9E3779B97F4A7C15ULL;
  for (int i = 0; i < 8; i++) {
    next_uniform(&state); // decorrelate small seeds
  }

  for (int p = 0; p < count; p++) {
    // radial: inverse CDF with linear interpolation inside the bin
    double u = next_uniform(&state) * total;
    int lo = 0;
    int hi = RADIAL_TABLE - 1;
    while (lo < hi) {
      int mid = (lo + hi) / 2;

      if (cdf[mid] < u) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }

    double prev = lo > 0 ? cdf[lo - 1] : 0.0;
    double frac = (cdf[lo] > prev) ? (u - prev) / (cdf[lo] - prev) : 0.5;
    double r = (lo + frac) * dr;

    // angular: rejection in (\cos(\theta, \phi))
    double theta = 0.0;
    double phi = 0.0;
    complex_t y = c_real(0.0);
    for (int tries = 0; tries < MAX_REJECT; tries++) {
      theta = acos(2.0 * next_uniform(&state) - 1.0);
      phi = 2.0 * M_PI * next_uniform(&state);
      y = harmonic(l, m, real_form, theta, phi);

      if (next_uniform(&state) * bound <= c_abs2(y)) {
        break;
      }
    }

    double st = sin(theta);
    double *o = out + 4 * (size_t)p;

    o[0] = r * st * cos(phi);
    o[1] = r * st * sin(phi);
    o[2] = r * cos(theta);

    double phase;
    if (real_form) {
      phase = (y.re * sign[lo] < 0.0) ? M_PI : 0.0;
    } else {
      phase = atan2(y.im, y.re);
      if (sign[lo] < 0) {
        phase += (phase > 0.0) ? -M_PI : M_PI;
      }
    }

    o[3] = phase;
  }

  free(cdf);
  free(sign);

  return 0;
}
