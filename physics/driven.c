/*
Driven two-level systems: general time-dependent Hamiltonians, Landau-Zener
sweeps, and lab-frame driving beyond the rotating-wave approximation (RWA).
*/

#include "driven.h"
#include "../core/complex.h"
#include "../core/vector.h"
#include <math.h>
#include <stdlib.h>

double time_fn_constant(double t, void *params) {
  return params ? *(double *)params : 0.0;
}

double time_fn_linear_ramp(double t, void *params) {
  double alpha = params ? *(double *)params : 0.0;

  return alpha * t;
}

// -i * H(t) * \psi for H(t) = (1/2)*[[\Delta(t), \Omega(t)],[\Omega(t),
// -\Delta(t)]]
static void rotating_frame_rhs(double t, const complex_t psi_in[2],
                               time_fn Delta, void *delta_params, time_fn Omega,
                               void *omega_params, complex_t psi_out[2]) {
  double D = Delta(t, delta_params);
  double O = Omega(t, omega_params);

  complex_t H_psi0 =
      c_add(c_scale(psi_in[0], 0.5 * D), c_scale(psi_in[1], 0.5 * O));
  complex_t H_psi1 =
      c_add(c_scale(psi_in[0], 0.5 * O), c_scale(psi_in[1], -0.5 * D));

  psi_out[0] = c_mul(c_imag(-1.0), H_psi0);
  psi_out[1] = c_mul(c_imag(-1.0), H_psi1);
}

int driven_two_level_evolve(cvector_t *psi, time_fn Delta, void *delta_params,
                            time_fn Omega, void *omega_params, double t0,
                            double dt, int steps) {
  if (!psi || psi->n != 2 || !Delta || !Omega || steps < 0) {
    return -1;
  }

  complex_t y[2] = {psi->data[0], psi->data[1]};
  double t = t0;

  for (int s = 0; s < steps; s++) {
    complex_t k1[2], k2[2], k3[2], k4[2], tmp[2];

    rotating_frame_rhs(t, y, Delta, delta_params, Omega, omega_params, k1);

    tmp[0] = c_add(y[0], c_scale(k1[0], dt * 0.5));
    tmp[1] = c_add(y[1], c_scale(k1[1], dt * 0.5));
    rotating_frame_rhs(t + dt * 0.5, tmp, Delta, delta_params, Omega,
                       omega_params, k2);

    tmp[0] = c_add(y[0], c_scale(k2[0], dt * 0.5));
    tmp[1] = c_add(y[1], c_scale(k2[1], dt * 0.5));
    rotating_frame_rhs(t + dt * 0.5, tmp, Delta, delta_params, Omega,
                       omega_params, k3);

    tmp[0] = c_add(y[0], c_scale(k3[0], dt));
    tmp[1] = c_add(y[1], c_scale(k3[1], dt));
    rotating_frame_rhs(t + dt, tmp, Delta, delta_params, Omega, omega_params,
                       k4);

    for (int i = 0; i < 2; i++) {
      complex_t sum = c_add(k1[i], c_scale(c_add(k2[i], k3[i]), 2.0));
      sum = c_add(sum, k4[i]);
      y[i] = c_add(y[i], c_scale(sum, dt / 6.0));
    }

    t += dt;
  }

  psi->data[0] = y[0];
  psi->data[1] = y[1];

  return 0;
}

double landau_zener_probability(double Omega, double alpha) {
  if (alpha <= 0.0) {
    return 0.0;
  }

  return exp(-M_PI * Omega * Omega / (2.0 * alpha));
}

// -i*H(t)*psi for the full lab-frame (non-RWA) Hamiltonian
static void lab_frame_rhs(double t, const complex_t psi_in[2], double omega0,
                          double Omega0, double omega_L, double phase,
                          complex_t psi_out[2]) {
  double drive = Omega0 * cos(omega_L * t + phase);

  complex_t H_psi0 =
      c_add(c_scale(psi_in[0], 0.5 * omega0), c_scale(psi_in[1], drive));
  complex_t H_psi1 =
      c_add(c_scale(psi_in[0], drive), c_scale(psi_in[1], -0.5 * omega0));

  psi_out[0] = c_mul(c_imag(-1.0), H_psi0);
  psi_out[1] = c_mul(c_imag(-1.0), H_psi1);
}

int driven_two_level_evolve_lab_frame(cvector_t *psi, double omega0,
                                      double Omega0, double omega_L,
                                      double phase, double t0, double dt,
                                      int steps) {
  if (!psi || psi->n != 2 || steps < 0) {
    return -1;
  }

  complex_t y[2] = {psi->data[0], psi->data[1]};
  double t = t0;

  for (int s = 0; s < steps; s++) {
    complex_t k1[2], k2[2], k3[2], k4[2], tmp[2];

    lab_frame_rhs(t, y, omega0, Omega0, omega_L, phase, k1);

    tmp[0] = c_add(y[0], c_scale(k1[0], dt * 0.5));
    tmp[1] = c_add(y[1], c_scale(k1[1], dt * 0.5));
    lab_frame_rhs(t + dt * 0.5, tmp, omega0, Omega0, omega_L, phase, k2);

    tmp[0] = c_add(y[0], c_scale(k2[0], dt * 0.5));
    tmp[1] = c_add(y[1], c_scale(k2[1], dt * 0.5));
    lab_frame_rhs(t + dt * 0.5, tmp, omega0, Omega0, omega_L, phase, k3);

    tmp[0] = c_add(y[0], c_scale(k3[0], dt));
    tmp[1] = c_add(y[1], c_scale(k3[1], dt));
    lab_frame_rhs(t + dt, tmp, omega0, Omega0, omega_L, phase, k4);

    for (int i = 0; i < 2; i++) {
      complex_t sum = c_add(k1[i], c_scale(c_add(k2[i], k3[i]), 2.0));
      sum = c_add(sum, k4[i]);
      y[i] = c_add(y[i], c_scale(sum, dt / 6.0));
    }

    t += dt;
  }

  psi->data[0] = y[0];
  psi->data[1] = y[1];

  return 0;
}
