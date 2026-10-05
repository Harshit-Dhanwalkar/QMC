/*
Time evolution of the 1D Dirac equation: split-operator Fourier method
*/

#include "dirac_evolve.h"
#include "../core/complex.h"
#include "../core/fft/fft.h"
#include "../core/vector.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int is_pow2_len(int n) { return n >= 4 && (n & (n - 1)) == 0; }

static int valid_spinor(const cvector_t *upper, const cvector_t *lower) {
  return upper && lower && upper->data && lower->data && upper->n == lower->n &&
         is_pow2_len(upper->n);
}

static int valid_constants(double dx, double m, double hbar, double c) {
  return isfinite(dx) && dx > 0.0 && isfinite(m) && m >= 0.0 && isfinite(hbar) &&
         hbar > 0.0 && isfinite(c) && c > 0.0;
}

// Wave number of momentum-space index j (same ordering as soft.c).
static double wave_number(int j, int n, double dx) {
  int jj = (j < n / 2) ? j : j - n;

  return 2.0 * M_PI * (double)jj / ((double)n * dx);
}

int dirac_evolve_1d(cvector_t *upper, cvector_t *lower, const double *V,
                    double dx, double dt, int steps, double m, double hbar,
                    double c) {
  if (!valid_spinor(upper, lower) || !valid_constants(dx, m, hbar, c) ||
      !isfinite(dt) || dt <= 0.0 || steps < 1) {
    return -1;
  }

  int n = upper->n;
  if (V) {
    for (int i = 0; i < n; i++) {
      if (!isfinite(V[i])) {
        return -1;
      }
    }
  }

  // Per-momentum quantities of the exact free propagator
  //   U(k) = cs - i sn (a \sigma_x + b \sigma_z),  a = c \hbar k,  b = m c^2
  // with cs = \cos(E dt/\hbar), sn = \sin(E dt/\hbar)/E (-> dt/\hbar as E -> 0)
  double *a = malloc((size_t)n * sizeof *a);
  double *cs = malloc((size_t)n * sizeof *cs);
  double *sn = malloc((size_t)n * sizeof *sn);
  complex_t *vphase = V ? malloc((size_t)n * sizeof *vphase) : NULL;
  cvector_t *work_u = cvector_alloc(n);
  cvector_t *work_l = cvector_alloc(n);
  if (!a || !cs || !sn || (V && !vphase) || !work_u || !work_l) {
    free(a);
    free(cs);
    free(sn);
    free(vphase);
    cvector_free(work_u);
    cvector_free(work_l);

    return -2;
  }

  double b = m * c * c;
  for (int j = 0; j < n; j++) {
    a[j] = c * hbar * wave_number(j, n, dx);
    double e = sqrt(a[j] * a[j] + b * b);
    double phase = e * dt / hbar;

    cs[j] = cos(phase);
    sn[j] = (e > 0.0) ? sin(phase) / e : dt / hbar;
  }

  if (V) {
    for (int i = 0; i < n; i++) {
      double phase = -V[i] * dt / (2.0 * hbar);

      vphase[i] = c_new(cos(phase), sin(phase));
    }
  }

  // Work on copies so a failure part-way cannot leave a half-evolved spinor
  for (int i = 0; i < n; i++) {
    work_u->data[i] = upper->data[i];
    work_l->data[i] = lower->data[i];
  }

  for (int s = 0; s < steps; s++) {
    if (V) {
      for (int i = 0; i < n; i++) {
        work_u->data[i] = c_mul(work_u->data[i], vphase[i]);
        work_l->data[i] = c_mul(work_l->data[i], vphase[i]);
      }
    }

    fft(work_u);
    fft(work_l);
    for (int j = 0; j < n; j++) {
      complex_t u = work_u->data[j];
      complex_t l = work_l->data[j];
      // (a sigma_x + b sigma_z)(u, l) = (b u + a l, a u - b l)
      complex_t hu = c_new(b * u.re + a[j] * l.re, b * u.im + a[j] * l.im);
      complex_t hl = c_new(a[j] * u.re - b * l.re, a[j] * u.im - b * l.im);

      // -i * sn * h = sn * (h.im, -h.re)
      work_u->data[j] =
          c_new(cs[j] * u.re + sn[j] * hu.im, cs[j] * u.im - sn[j] * hu.re);
      work_l->data[j] =
          c_new(cs[j] * l.re + sn[j] * hl.im, cs[j] * l.im - sn[j] * hl.re);
    }

    ifft(work_u);
    ifft(work_l);

    if (V) {
      for (int i = 0; i < n; i++) {
        work_u->data[i] = c_mul(work_u->data[i], vphase[i]);
        work_l->data[i] = c_mul(work_l->data[i], vphase[i]);
      }
    }
  }

  for (int i = 0; i < n; i++) {
    upper->data[i] = work_u->data[i];
    lower->data[i] = work_l->data[i];
  }

  free(a);
  free(cs);
  free(sn);
  free(vphase);
  cvector_free(work_u);
  cvector_free(work_l);

  return 0;
}

int dirac_packet_1d(cvector_t *upper, cvector_t *lower, double dx, double x0,
                    double k0, double sigma, int branch, double m, double hbar,
                    double c) {
  if (!valid_spinor(upper, lower) || !valid_constants(dx, m, hbar, c) ||
      !isfinite(x0) || !isfinite(k0) || !isfinite(sigma) || sigma <= 0.0 ||
      branch < -1 || branch > 1) {
    return -1;
  }

  int n = upper->n;
  double b = m * c * c;
  double a0 = c * hbar * k0;
  double e0 = sqrt(a0 * a0 + b * b);

  // Spinor of the requested branch at k0: eigenvector (E_s + b, a0) of H_0(k0)
  // with eigenvalue E_s = branch * E_k0 (any fixed spinor for branch 0)
  double s_up = 1.0;
  double s_lo = 0.0;
  if (branch != 0) {
    double es = (double)branch * e0;
    double norm = hypot(es + b, a0);

    if (norm > 1e-12 * (e0 + b + 1.0)) {
      s_up = (es + b) / norm;
      s_lo = a0 / norm;
    } else if (branch < 0) {
      s_up = 0.0; // k0 = 0: negative-energy state is the lower component
      s_lo = 1.0;
    }
  }

  for (int i = 0; i < n; i++) {
    double x = ((double)i - 0.5 * (double)n) * dx - x0;
    double env = exp(-x * x / (4.0 * sigma * sigma));
    complex_t g = c_new(env * cos(k0 * x), env * sin(k0 * x));

    upper->data[i] = c_scale(g, s_up);
    lower->data[i] = c_scale(g, s_lo);
  }

  if (branch != 0) {
    // Project onto the sign-`branch` energy subspace in momentum space:
    // P = (1 + branch H_0(k)/E_k) / 2. Needed because a Gaussian of finite
    // width mixes in neighbouring k whose eigenspinors differ from the one
    // at k0
    fft(upper);
    fft(lower);
    for (int j = 0; j < n; j++) {
      double a = c * hbar * wave_number(j, n, dx);
      double e = sqrt(a * a + b * b);
      complex_t u = upper->data[j];
      complex_t l = lower->data[j];

      if (e > 0.0) {
        double f = (double)branch / e;
        complex_t hu = c_new(b * u.re + a * l.re, b * u.im + a * l.im);
        complex_t hl = c_new(a * u.re - b * l.re, a * u.im - b * l.im);

        upper->data[j] = c_new(0.5 * (u.re + f * hu.re), 0.5 * (u.im + f * hu.im));
        lower->data[j] = c_new(0.5 * (l.re + f * hl.re), 0.5 * (l.im + f * hl.im));
      }
    }

    ifft(upper);
    ifft(lower);
  }

  double norm = dirac_norm_1d(upper, lower, dx);
  if (!(norm > 0.0) || !isfinite(norm)) {
    return -1;
  }

  double scale = 1.0 / sqrt(norm);
  for (int i = 0; i < n; i++) {
    upper->data[i] = c_scale(upper->data[i], scale);
    lower->data[i] = c_scale(lower->data[i], scale);
  }

  return 0;
}

double dirac_norm_1d(const cvector_t *upper, const cvector_t *lower,
                     double dx) {
  if (!upper || !lower || !upper->data || !lower->data ||
      upper->n != lower->n || upper->n < 1 || !isfinite(dx) || dx <= 0.0) {
    return NAN;
  }

  double sum = 0.0;
  for (int i = 0; i < upper->n; i++) {
    sum += c_abs2(upper->data[i]) + c_abs2(lower->data[i]);
  }

  return sum * dx;
}

double dirac_position_1d(const cvector_t *upper, const cvector_t *lower,
                         double dx) {
  double norm = dirac_norm_1d(upper, lower, dx);
  if (!(norm > 0.0)) {
    return NAN;
  }

  int n = upper->n;
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    double x = ((double)i - 0.5 * (double)n) * dx;

    sum += x * (c_abs2(upper->data[i]) + c_abs2(lower->data[i]));
  }

  return sum * dx / norm;
}

double dirac_energy_1d(const cvector_t *upper, const cvector_t *lower,
                       const double *V, double dx, double m, double hbar,
                       double c) {
  double norm = dirac_norm_1d(upper, lower, dx);
  if (!(norm > 0.0) || !is_pow2_len(upper->n) || !valid_constants(dx, m, hbar, c)) {
    return NAN;
  }

  int n = upper->n;
  cvector_t *u = cvector_copy(upper);
  cvector_t *l = cvector_copy(lower);
  if (!u || !l) {
    cvector_free(u);
    cvector_free(l);

    return NAN;
  }

  fft(u);
  fft(l);

  double b = m * c * c;
  double kinetic = 0.0;
  for (int j = 0; j < n; j++) {
    double a = c * hbar * wave_number(j, n, dx);
    complex_t uj = u->data[j];
    complex_t lj = l->data[j];

    // Re[ u* (b u + a l) + l* (a u - b l) ]
    kinetic += b * (c_abs2(uj) - c_abs2(lj)) +
               2.0 * a * (uj.re * lj.re + uj.im * lj.im);
  }

  // Parseval for unnormalised forward FFT: ]sum_i |f_i|^2 = (1/N) \sum_j |F_j|^2
  kinetic *= dx / (double)n;

  double potential = 0.0;
  if (V) {
    for (int i = 0; i < n; i++) {
      potential += V[i] * (c_abs2(upper->data[i]) + c_abs2(lower->data[i]));
    }

    potential *= dx;
  }

  cvector_free(u);
  cvector_free(l);

  return (kinetic + potential) / norm;
}

double dirac_step_transmission(double E, double V0, double m, double hbar,
                               double c) {
  if (!isfinite(E) || !isfinite(V0) || !isfinite(m) || m < 0.0 ||
      !isfinite(hbar) || hbar <= 0.0 || !isfinite(c) || c <= 0.0) {
    return NAN;
  }

  double b = m * c * c;
  if (E <= b) {
    return NAN;
  }

  double ep = E - V0;
  if (fabs(ep) <= b) {
    return 0.0; // no propagating state inside the step
  }

  double k = sqrt(E * E - b * b) / (hbar * c);
  double q = copysign(sqrt(ep * ep - b * b) / (hbar * c), ep);
  double kappa = (q / k) * (E + b) / (ep + b);

  return 4.0 * kappa / ((1.0 + kappa) * (1.0 + kappa));
}
