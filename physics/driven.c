/*
Driven two-level systems: general time-dependent Hamiltonians, Landau-Zener
sweeps, and lab-frame driving beyond the rotating-wave approximation (RWA).
*/

#include "driven.h"
#include "../core/complex.h"
#include "../core/ode/rk4.h"
#include "../core/vector.h"
#include <math.h>
#include <stdlib.h>

static const double HALF = 0.5;

double time_fn_constant(double time, void *params) {
  return params ? *(double *)params : 0.0;
}

double time_fn_linear_ramp(double time, void *params) {
  double alpha = params ? *(double *)params : 0.0;

  return alpha * time;
}

// Context for rotating-frame RHS
typedef struct {
  time_fn Delta;
  void *delta_params;
  time_fn Omega;
  void *omega_params;
} driven_rwa_ctx_t;

// -i * H(t) * \psi for
// H(t) = (1/2) * [[\Delta(t), \Omega(t)],[\Omega(t), -\Delta(t)]]
static void rotating_frame_rhs(double time, const cvector_t *state,
                               cvector_t *dydt, void *params) {
  driven_rwa_ctx_t *ctx = (driven_rwa_ctx_t *)params;
  double detuning = ctx->Delta(time, ctx->delta_params);
  double rabi_freq = ctx->Omega(time, ctx->omega_params);

  complex_t H_psi0 = c_add(c_scale(state->data[0], HALF * detuning),
                           c_scale(state->data[1], HALF * rabi_freq));
  complex_t H_psi1 = c_add(c_scale(state->data[0], HALF * rabi_freq),
                           c_scale(state->data[1], -HALF * detuning));

  dydt->data[0] = c_mul(c_imag(-1.0), H_psi0);
  dydt->data[1] = c_mul(c_imag(-1.0), H_psi1);
}

int driven_two_level_evolve(cvector_t *psi, time_fn Delta, void *delta_params,
                            time_fn Omega, void *omega_params,
                            driven_params_t params) {
  if (!psi || psi->n != 2 || !Delta || !Omega || params.steps < 0) {
    return -1;
  }

  driven_rwa_ctx_t ctx = {Delta, delta_params, Omega, omega_params};
  double time = params.t0;

  for (int s = 0; s < params.steps; s++) {
    if (rk4_step(time, params.dt, psi, rotating_frame_rhs, &ctx) != 0) {
      return -1;
    }

    time += params.dt;
  }

  return 0;
}

double landau_zener_probability(double Omega, double alpha) {
  if (alpha <= 0.0) {
    return 0.0;
  }

  return exp(-M_PI * Omega * Omega / (2.0 * alpha));
}

// Context for lab-frame (non-RWA) RHS.
typedef struct {
  double omega0;
  double Omega0;
  double omega_L;
  double phase;
} driven_lab_ctx_t;

// -i * H(t) * \psi for full lab-frame (non-RWA) Hamiltonian
static void lab_frame_rhs(double t, const cvector_t *y, cvector_t *dydt,
                          void *params) {
  const driven_lab_ctx_t *c = (driven_lab_ctx_t *)params;
  double drive = c->Omega0 * cos(c->omega_L * t + c->phase);

  complex_t H_psi0 =
      c_add(c_scale(y->data[0], 0.5 * c->omega0), c_scale(y->data[1], drive));
  complex_t H_psi1 =
      c_add(c_scale(y->data[0], drive), c_scale(y->data[1], -0.5 * c->omega0));

  dydt->data[0] = c_mul(c_imag(-1.0), H_psi0);
  dydt->data[1] = c_mul(c_imag(-1.0), H_psi1);
}

int driven_two_level_evolve_lab_frame(cvector_t *psi, double omega0,
                                      double Omega0, double omega_L,
                                      double phase, double t0, double dt,
                                      int steps) {
  if (!psi || psi->n != 2 || steps < 0) {
    return -1;
  }

  driven_lab_ctx_t ctx = {omega0, Omega0, omega_L, phase};
  double time = t0;

  for (int s = 0; s < steps; s++) {
    if (rk4_step(time, dt, psi, lab_frame_rhs, &ctx) != 0) {
      return -1;
    }

    time += dt;
  }

  return 0;
}
