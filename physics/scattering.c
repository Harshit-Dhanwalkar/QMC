/*
Scattering: phase shifts and Born approximation.
 TODO: cross-section
*/

#include "scattering.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/ode/numerov.h"
#include "../core/special/special.h"
#include "../core/vector.h"
#include "potentials.h"
#include <math.h>
#include <stdlib.h>

// Spherical Bessel/Neumann functions via upward recurrence
static double spherical_bessel_j(int l, double x) {
  if (l == 0)
    return sin(x) / x;
  double j0 = sin(x) / x;
  double j1 = sin(x) / (x * x) - cos(x) / x;
  if (l == 1)
    return j1;
  double jm1 = j0, j = j1;
  for (int n = 1; n < l; n++) {
    double jp1 = (2 * n + 1) / x * j - jm1;
    jm1 = j;
    j = jp1;
  }
  return j;
}

static double spherical_neumann_y(int l, double x) {
  if (l == 0)
    return -cos(x) / x;
  double y0 = -cos(x) / x;
  double y1 = -cos(x) / (x * x) - sin(x) / x;
  if (l == 1)
    return y1;
  double ym1 = y0, yv = y1;
  for (int n = 1; n < l; n++) {
    double yp1 = (2 * n + 1) / x * yv - ym1;
    ym1 = yv;
    yv = yp1;
  }
  return yv;
}

// phase_shift: scattering phase shift for partial wave l via integration +
// asymptotic matching.
/*
 * Method: integrate radial equation for
 * u_l(r) = r * R_l(r),
 * -hbar_sq_2m u_l'' + [V(r) + hbar_sq_2m * l(l+1) / r^2] u_l = E u_l,
 * E = hbar_sq_2m * k^2, using
 *
 *   u(r) ~ C*[\cos(\delta) * (kr * j_l(kr)) - \sin(\delta) * (kr * n_l(kr))]
 *
 * Solving resulting 2x2 system for \delta :
 *   \tan(\delta) = (u2*J1 - u1*J2) / (u2*N1 - u1*N2)
 * where J,N are Riccati-Bessel/Neumann functions kr * j_l(kr), kr * n_l(kr)
 * evaluated at r1,r2.
 */
double phase_shift(int l, double k, potential_fn V, void *params, double r_min,
                   double r_max, int N, double hbar_sq_2m) {
  if (!V || N < 100 || k <= 0.0 || r_min <= 0.0 || r_max <= r_min)
    return 0.0;

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
  double J1 = kr1 * spherical_bessel_j(l, kr1);
  double Nn1 = kr1 * spherical_neumann_y(l, kr1);
  double J2 = kr2 * spherical_bessel_j(l, kr2);
  double Nn2 = kr2 * spherical_neumann_y(l, kr2);

  double num = u2 * J1 - u1 * J2;
  double den = u2 * Nn1 - u1 * Nn2;
  return atan2(num, den);
}

complex_t born_amplitude(potential_fn V, void *params, double k, double theta,
                         double r_max, int N) {
  // Born approximation: f(θ) = - (2m)/(4\pi \hbar^2) ∫ V(r) e^{-i q·r} d³r
  // For central potential: f(θ) = - (2m/\hbar^2) (1/q) ∫_0^∞ r V(r) sin(q r) dr
  // where q = 2k sin(θ/2).
  // Integrate numerically
  if (!V || N < 2)
    return c_zero();
  double q = 2.0 * k * sin(theta / 2.0);
  double dr = r_max / (N - 1);
  double integral = 0.0;
  for (int i = 0; i < N; i++) {
    double r = i * dr;
    if (r < 1e-15)
      continue;
    double Vr = V(r, params);
    integral += r * Vr * sin(q * r) * dr;
  }
  complex_t f = c_real(-2.0 * M_ELECTRON / (HBAR * HBAR) * integral / q);
  return f;
}

double born_cross_section(complex_t f_theta) { return c_abs2(f_theta); }
