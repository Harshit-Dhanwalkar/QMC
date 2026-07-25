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
