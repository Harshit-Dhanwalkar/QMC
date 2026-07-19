/*
Two-level Rabi oscillations (rotating-wave approximation).
*/

#include "rabi.h"
#include "../core/complex.h"
#include "../core/ode/crank_nicolson.h"
#include "../core/vector.h"
#include <math.h>

double rabi_excited_probability(double t, double Omega, double Delta) {
  double Omega_R = sqrt(Omega * Omega + Delta * Delta);
  if (Omega_R < 1e-300)
    return 0.0; // no coupling and no detuning: nothing happens
  double s = sin(0.5 * Omega_R * t);

  return (Omega * Omega / (Omega_R * Omega_R)) * s * s;
}

int rabi_evolve_exact(cvector_t *psi, double t, double Omega, double Delta) {
  if (!psi || psi->n != 2)
    return -1;

  double Omega_R = sqrt(Omega * Omega + Delta * Delta);
  double half = 0.5 * Omega_R * t;
  double cos_h = cos(half);
  // \sinc_h = \sin(half) / \Omega_R, with the \Omega_R->0 limit (\sin(x)/x->1)
  double sinc_h = (Omega_R > 1e-300) ? sin(half) / Omega_R : 0.5 * t;

  // U = \cos_h*I - i * \sinc_h * K,
  // K = [[\Delta, \Omega],[\Omega, -\Delta]]
  complex_t U11 = c_add(c_real(cos_h), c_imag(-sinc_h * Delta));
  complex_t U12 = c_imag(-sinc_h * Omega);
  complex_t U21 = c_imag(-sinc_h * Omega);
  complex_t U22 = c_add(c_real(cos_h), c_imag(sinc_h * Delta));

  complex_t p0 = psi->data[0], p1 = psi->data[1];
  psi->data[0] = c_add(c_mul(U11, p0), c_mul(U12, p1));
  psi->data[1] = c_add(c_mul(U21, p0), c_mul(U22, p1));

  return 0;
}

int rabi_evolve_numerical(cvector_t *psi, double hbar, double Omega,
                          double Delta, double dt, int steps) {
  if (!psi || psi->n != 2 || steps < 1)
    return -1;

  double diag[2] = {0.5 * hbar * Delta, -0.5 * hbar * Delta};
  double offdiag[1] = {0.5 * hbar * Omega};

  for (int s = 0; s < steps; s++) {
    if (crank_nicolson_step(diag, offdiag, dt, psi) != 0)
      return -1;
  }

  return 0;
}
