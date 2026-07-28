/*
Scattering: phase shifts, Born approximation, partial-wave observables,
and Lippmann-Schwinger / Born-series extension.
*/

#include "scattering.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/linalg/lu.h"
#include "../core/matrix.h"
#include "../core/ode/numerov.h"
#include "../core/special/special.h"
#include "../core/vector.h"
#include "potentials.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// phase_shift: scattering phase shift for partial wave l via integration +
// asymptotic matching.
/*
 * Method: integrate radial equation for
 * u_l(r) = r * R_l(r),
 * -hbar_sq_2m u_l'' + [V(r) + hbar_sq_2m * l(l+1) / r^2] u_l = E u_l,
 * E = hbar_sq_2m * k^2, using
 *
 *   u(r) ~ C*[\cos(\delta) * S_l(kr) - \sin(\delta) * C_l(kr)]
 * Where S_l, C_l are the Riccati-Bessel functions evaluated at r1, r2.
 *
 * Solving resulting 2x2 system for \delta :
 *   \tan(\delta) = (u2*J1 - u1*J2) / (u2 * N1 - u1 * N2)
 * Where J,N are Riccati-Bessel functions S_l, C_l evaluated at r1,r2. j_l, n_l
 * come from core/special (sph_bessel_j/ sph_bessel_y), which use a numerically
 * stable downward (Miller's algorithm) recurrence for j_l
 */
double phase_shift(int l, double k, potential_fn V, void *params, double r_min,
                   double r_max, int N, double hbar_sq_2m) {
  if (!V || N < 100 || k <= 0.0 || r_min <= 0.0 || r_max <= r_min) {
    return 0.0;
  }

  double dr = (r_max - r_min) / (N - 1);
  double *r = malloc(N * sizeof *r);
  double *V_eff = malloc(N * sizeof *V_eff);
  if (!r || !V_eff) {
    free(r);
    free(V_eff);

    return 0.0;
  }

  for (int i = 0; i < N; i++) {
    r[i] = r_min + i * dr;
    V_eff[i] = V(r[i], params) + hbar_sq_2m * l * (l + 1.0) / (r[i] * r[i]);
  }

  numerov_params_t np = {
      .x = r, .V = V_eff, .n = N, .dx = dr, .hbar_sq_2m = hbar_sq_2m};
  cvector_t *u = cvector_alloc(N);
  if (!u) {
    free(r);
    free(V_eff);

    return 0.0;
  }

  double E = hbar_sq_2m * k * k;
  numerov_integrate(&np, E, u);

  // 2 matching points well into asymptotic region
  int i1 = N - 100 > 0 ? N - 100 : N / 2;
  int i2 = N - 10;
  double r1 = r[i1], r2 = r[i2];
  double u1 = u->data[i1].re, u2 = u->data[i2].re;

  free(r);
  free(V_eff);
  cvector_free(u);

  double kr1 = k * r1, kr2 = k * r2;
  double J1 = riccati_bessel_j(l, kr1);
  double Nn1 = riccati_bessel_y(l, kr1);
  double J2 = riccati_bessel_j(l, kr2);
  double Nn2 = riccati_bessel_y(l, kr2);

  double num = u2 * J1 - u1 * J2;
  double den = u2 * Nn1 - u1 * Nn2;
  return atan2(num, den);
}

complex_t born_amplitude(potential_fn V, void *params, double k, double theta,
                         double r_max, int N, double hbar_sq_2m) {
  // Born approximation: f(\theta) = - (2m)/(4 *\pi * \hbar^2) \int V(r) * e^{-i
  // q \cdot r} d^3r
  // For central potential: f(\theta) = - (2m/\hbar^2) (1/q)
  // \int_0^\intfy r * V(r) * \sin(q*r) dr where q = 2k * \sin(\theta/2).
  // Integrate numerically
  if (!V || N < 2 || hbar_sq_2m <= 0.0) {
    return c_zero();
  }

  double q = 2.0 * k * sin(theta / 2.0);
  double dr = r_max / (N - 1);
  double integral = 0.0;
  for (int i = 0; i < N; i++) {
    double r = i * dr;
    if (r < 1e-15) {
      continue;
    }

    double Vr = V(r, params);
    integral += r * Vr * sin(q * r) * dr;
  }
  complex_t f = c_real(-integral / (hbar_sq_2m * q));

  return f;
}

double born_cross_section(complex_t f_theta) { return c_abs2(f_theta); }

/* Partial-wave observables */
complex_t s_matrix_element(double delta_l) {
  return c_new(cos(2.0 * delta_l), sin(2.0 * delta_l));
}

complex_t partial_wave_amplitude(const double *delta_l, int l_max, double k,
                                 double theta) {
  if (!delta_l || l_max < 0 || k <= 0.0) {
    return c_zero();
  }

  double x = cos(theta);
  complex_t f = c_zero();
  for (int l = 0; l <= l_max; l++) {
    double Pl = legendre(l, x);
    double sd = sin(delta_l[l]);
    complex_t phase = c_new(cos(delta_l[l]), sin(delta_l[l]));
    complex_t term = c_scale(phase, (2.0 * l + 1.0) * sd * Pl);
    f = c_add(f, term);
  }

  return c_scale(f, 1.0 / k);
}

double total_cross_section_partial_waves(const double *delta_l, int l_max,
                                         double k) {
  if (!delta_l || l_max < 0 || k <= 0.0) {
    return 0.0;
  }

  double sigma = 0.0;
  for (int l = 0; l <= l_max; l++) {
    double sd = sin(delta_l[l]);
    sigma += (2.0 * l + 1.0) * sd * sd;
  }

  return 4.0 * M_PI / (k * k) * sigma;
}

/* Lippmann-Schwinger / Born series */
// Radial LS equation
// Riccati-Neumann sign convention (riccati_bessel_y = x*sph_bessel_y(x)
//   u_l(r) = hatj_l(kr) + coeff * \int_0^\infty hatj_l(kr_<) hatn_l(kr_>) V(r')
//   u_l(r') dr'
//   \tan(\delta_l) = -coeff * \int_0^\infty hatj_l(kr) V(r) u_l(r) dr
//
// Where coeff = 1/(hbar_sq_2m * k). Setting u_l = hatj_l(kr) (for no
// self-consistency) reduces \tan(\delta_l) formula to first-Born partial-wave
// phase shift.

static int build_free_wave_grid(double r_min, double r_max, int N, double *r,
                                double *w) {
  if (N < 2) {
    return -1;
  }

  double dr = (r_max - r_min) / (N - 1);
  for (int i = 0; i < N; i++) {
    r[i] = r_min + i * dr;
    w[i] = dr;
  }
  w[0] *= 0.5;
  w[N - 1] *= 0.5;

  return 0;
}

double phase_shift_lippmann_schwinger(int l, double k, potential_fn V,
                                      void *params, double r_min, double r_max,
                                      int N, double hbar_sq_2m) {
  if (!V || N < 2 || k <= 0.0 || r_min < 0.0 || r_max <= r_min ||
      hbar_sq_2m <= 0.0) {
    return 0.0;
  }

  double *r = malloc(N * sizeof *r);
  double *w = malloc(N * sizeof *w);
  double *Vr = malloc(N * sizeof *Vr);
  double *jkr = malloc(N * sizeof *jkr);
  if (!r || !w || !Vr || !jkr || build_free_wave_grid(r_min, r_max, N, r, w)) {
    free(r);
    free(w);
    free(Vr);
    free(jkr);

    return 0.0;
  }

  for (int i = 0; i < N; i++) {
    Vr[i] = V(r[i], params);
    jkr[i] = riccati_bessel_j(l, k * r[i]);
  }

  cmatrix_t *M = cmatrix_alloc(N, N);
  cvector_t *b = cvector_alloc(N);
  cvector_t *u = cvector_alloc(N);
  if (!M || !b || !u) {
    free(r);
    free(w);
    free(Vr);
    free(jkr);
    if (M)
      cmatrix_free(M);
    if (b)
      cvector_free(b);
    if (u)
      cvector_free(u);

    return 0.0;
  }

  double coeff = 1.0 / (hbar_sq_2m * k);
  for (int i = 0; i < N; i++) {
    b->data[i] = c_real(jkr[i]);
    for (int j = 0; j < N; j++) {
      double rlo = r[i] < r[j] ? r[i] : r[j];
      double rhi = r[i] < r[j] ? r[j] : r[i];
      double G = riccati_bessel_j(l, k * rlo) * riccati_bessel_y(l, k * rhi);
      double Aij = coeff * w[j] * G * Vr[j];
      double diag = (i == j) ? 1.0 : 0.0;
      CMAT(M, i, j) = c_real(diag - Aij);
    }
  }

  double tan_delta = 0.0;
  int *pivot = lu_decompose(M);
  if (pivot) {
    if (lu_solve(M, pivot, b, u) == 0) {
      double integral = 0.0;
      for (int i = 0; i < N; i++) {
        integral += w[i] * jkr[i] * Vr[i] * u->data[i].re;
      }
      tan_delta = -coeff * integral;
    }
    free(pivot);
  }

  cmatrix_free(M);
  cvector_free(b);
  cvector_free(u);
  free(r);
  free(w);
  free(Vr);
  free(jkr);

  return atan(tan_delta);
}

double born_series_phase_shift(int l, double k, potential_fn V, void *params,
                               double r_min, double r_max, int N,
                               double hbar_sq_2m, int order) {
  if (!V || N < 2 || k <= 0.0 || r_min < 0.0 || r_max <= r_min ||
      hbar_sq_2m <= 0.0 || order < 1) {
    return 0.0;
  }

  double *r = malloc(N * sizeof *r);
  double *w = malloc(N * sizeof *w);
  double *Vr = malloc(N * sizeof *Vr);
  double *jkr = malloc(N * sizeof *jkr);
  double *u = malloc(N * sizeof *u);
  double *u_next = malloc(N * sizeof *u_next);
  if (!r || !w || !Vr || !jkr || !u || !u_next ||
      build_free_wave_grid(r_min, r_max, N, r, w)) {
    free(r);
    free(w);
    free(Vr);
    free(jkr);
    free(u);
    free(u_next);

    return 0.0;
  }

  for (int i = 0; i < N; i++) {
    Vr[i] = V(r[i], params);
    jkr[i] = riccati_bessel_j(l, k * r[i]);
    u[i] = jkr[i]; // zeroth Neumann iterate: u^(0) = hatj_l(kr)
  }

  double coeff = 1.0 / (hbar_sq_2m * k);

  // order-1 iterations of the LS kernel
  for (int n = 1; n < order; n++) {
    for (int i = 0; i < N; i++) {
      double acc = 0.0;
      for (int j = 0; j < N; j++) {
        double rlo = r[i] < r[j] ? r[i] : r[j];
        double rhi = r[i] < r[j] ? r[j] : r[i];
        double G = riccati_bessel_j(l, k * rlo) * riccati_bessel_y(l, k * rhi);
        acc += w[j] * G * Vr[j] * u[j];
      }
      u_next[i] = jkr[i] + coeff * acc;
    }
    memcpy(u, u_next, N * sizeof *u);
  }

  double integral = 0.0;
  for (int i = 0; i < N; i++) {
    integral += w[i] * jkr[i] * Vr[i] * u[i];
  }
  double tan_delta = -coeff * integral;

  free(r);
  free(w);
  free(Vr);
  free(jkr);
  free(u);
  free(u_next);

  return atan(tan_delta);
}
