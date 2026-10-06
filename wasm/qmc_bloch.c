/*
 * qmc_bloch.c - WebAssembly front end for Bloch-sphere demo
 *
 * Physics: physics/bloch.c (RK4 Lindblad evolution of a driven, damped qubit,
 * physics/lindblad.c underneath) and physics/rabi.c (closed-form Rabi curve
 * for comparison). This file only keeps current state
 */

#include "../physics/bloch.h"
#include "../physics/rabi.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

#define BL_DT 0.01

static double g_v[3] = {0.0, 0.0, 1.0};
static double g_t;

QMC_EXPORT(qmc_bl_dt)
double qmc_bl_dt(void) { return BL_DT; }

/* Pointer to 3 doubles (x, y, z) of Bloch vector */
QMC_EXPORT(qmc_bl_vec)
double *qmc_bl_vec(void) { return g_v; }

QMC_EXPORT(qmc_bl_time)
double qmc_bl_time(void) { return g_t; }

/* Set state (|v| <= 1) and restart clock. Returns 0 on success */
QMC_EXPORT(qmc_bl_reset)
int qmc_bl_reset(double x, double y, double z) {
  if (!isfinite(x) || !isfinite(y) || !isfinite(z) ||
      sqrt(x * x + y * y + z * z) > 1.0 + 1e-9) {
    return -1;
  }
  g_v[0] = x;
  g_v[1] = y;
  g_v[2] = z;
  g_t = 0.0;

  return 0;
}

/* Advance `steps` steps of 0.01. Returns 0 on success */
QMC_EXPORT(qmc_bl_step)
int qmc_bl_step(double omega, double delta, double gamma1, double gamma_phi,
                int steps) {
  int rc = bloch_evolve(g_v, omega, delta, gamma1, gamma_phi, BL_DT, steps);

  if (rc == 0) {
    g_t += BL_DT * steps;
  }

  return rc;
}

/* Purity Tr(rho^2) = (1 + |v|^2)/2. */
QMC_EXPORT(qmc_bl_purity)
double qmc_bl_purity(void) {
  return 0.5 * (1.0 + g_v[0] * g_v[0] + g_v[1] * g_v[1] + g_v[2] * g_v[2]);
}

/* Closed-form excited-state probability from ground state. */
QMC_EXPORT(qmc_bl_rabi)
double qmc_bl_rabi(double t, double omega, double delta) {
  return rabi_excited_probability(t, omega, delta);
}
