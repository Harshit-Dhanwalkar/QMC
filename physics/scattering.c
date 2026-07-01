/*
Scattering: phase shifts and Born approximation.
 TODO: cross-section
*/
#include "scattering.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/ode/numerov.h"
#include "../core/special/special.h"
#include "potentials.h"
#include <math.h>
#include <stdlib.h>

double phase_shift(int l, double k, potential_fn V, void *params, double r_min,
                   double r_max, int N) {
  // Solve radial Schrödinger for given l and k (energy E = \hbar^2k^2/2m)
  // Compare asymptotic solution to spherical Bessel/Neumann.
  // HACK: Implement simple integration (Numerov) and extract phase shift.
  // return 0 as placeholder.
  return 0.0;
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
