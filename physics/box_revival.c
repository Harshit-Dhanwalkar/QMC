/*
 * Quantum revivals in a particle-in-a-box
 *
 * Units: Well is [0, 1], \hbar = 1 and mass is chosen so that eigenenergies are
 * E_n = n^2
 * Eigenfunctions are
 *     \phi_n(x) = \sqrt(2) \sin(n \pi x)
 * so time-evolved state is
 *     \psi(x, t) = \sum_{n >= 1} c_n \exp(-i n^2 t) \phi_n(x)
 * and whole system is periodic with revival time T = 2 \pi (harmonic-oscillator
 * analogy is misleading; "clock" here is n^2 phase, which returns to 1 for
 * every integer n only at multiples of 2 \pi)
 *
 * Projection (box_gaussian_coefficients):
 *   Initial packet is a Gaussian of spread sigma centred at x0,multiplied by a
 *   plane-wave phase e^{i k0 x} with k0 = \pi n0 (so n0 is central mode number)
 *   Overlap integral
 *       c_n = \int_0^1 \sqrt(2) \sin(n \pi x)
 *                      \exp^{-(x - x0)^2 / (4 \sigma^2)} \exp^{i k0 x} dx
 *   is evaluated by composite Simpson's rule on [0, 1]. Number of subintervals
 *   is scaled with highest mode number requested plus central mode, since
 *   integrand oscillates roughly (n + n0) half-waves over interval and Simpson
 *   grid has to resolve them
 *
 *   Gaussian is not normalised on [0, 1] (only on real line), so after
 *   coefficients are computed routine computes on-grid norm of unprojected
 *   packet and divides every coefficient through by its square root. This fixes
 *   amplitude to 1 exactly, independent of how much of Gaussian actually sat
 *   inside well. Retained weight (sum |c_n|^2 before that second
 *   renormalisation) is returned so caller can tell how much of packet fell
 *   outside truncated mode space; a value near 1 means truncation is not losing
 *   weight
 *
 * Time evolution (box_wavefunction, box_density_map, box_overlap):
 *   All three routines fold t modulo T = 2 \pi before touching phase and inside
 *   that they fold n^2 t modulo 2 \pi per mode. This is what keeps phase
 *   accurate for t well beyond T: without reduction, argument to exp() would
 *   grow as n^2 t and lose precision long before interesting revival structure
 *   appears
 *
 *   box_wavefunction returns \psi(x_j, t) at nx supplied positions, with
 *   eigenfunction phases precomputed into a stack array of at most
 *   BOX_MAX_MODES complex doubles
 *
 *   box_density_map fills a full space-time grid \rho[i * nx + j] for
 *   nt evenly spaced times from 0 to t_max. It is same computation as
 *   calling box_wavefunction nt times, but it precomputes sine table
 *   \sin(n \pi x_j) once (this is dominant cost), and it skips modes
 *   whose |c_n|^2 is below 1e-16, so a sharply peaked packet pays for only
 *   modes it actually occupies. Uses one heap allocation of (n_modes + 2) * nx
 *   doubles
 *
 *   box_overlap returns autocorrelation function
 *       <\psi(0)|\psi(t)> = \sum_n |c_n|^2 e^{-i n^2 t},
 *   which starts at 1 and passes through a sequence of partial revivals
 *   ("carpet") that mirrors Talbot effect in optics. Because phase prefactor is
 *   |c_n|^2 rather than c_n, this is return probability amplitude for a
 *   measurement in initial state and |<\psi(0)|\psi(t)>|^2
 *
 * Closed-form references (box_mean_energy, box_mean_position):
 *    mean energy  = \sum_n n^2 |c_n|^2
 *    mean position = \int x |\psi(x)|^2 dx / \int |\psi(x)|^2 dx
 *  second computed on caller's x grid by trapezoid rule on both numerator and
 *  denominator so ratio is well defined even for an unnormalised grid sample
 *
 * Edge cases:
 *   x0 must be strictly inside well (0 < x0 < 1) and sigma confined to
 *   [0.005, 0.25]
 *   outside that range projection either misses well or starts fighting grid
 *   resolution. n_modes is capped at BOX_MAX_MODES (512) so stack array in
 *   box_wavefunction stays bounded and coefficient buffer fits comfortably in
 *   cache
 */

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
  // Round subinterval count up to an even number so
  // Simpson's rule is exact for odd/even weighting pattern
  int steps = 8 * (n_modes + (int)fabs(n0) + 16);
  steps += steps % 2;

  double h = 1.0 / steps;
  double k0 = M_PI * n0;
  double norm = 0.0;

  for (int n = 1; n <= n_modes; n++) {
    double re = 0.0;
    double im = 0.0;

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
