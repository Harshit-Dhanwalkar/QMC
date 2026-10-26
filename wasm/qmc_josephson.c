/*
 * qmc_josephson.c - WebAssembly front end for Josephson-junction demo
 *
 * Physics: physics/josephson.c (RCSJ model: damped driven pendulum in a tilted
 * washboard potential). Time is in units of (1 / \omega_p) and voltage is
 * <\phi'>
 * A live state is stepped by page; I-V curves (up and down sweeps) are computed
 * for current Q, ac amplitude and frequency
 */

#include "../physics/josephson.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { JJ_CURVE = 161 };

static jj_params_t g_p = {2.0, 0.5, 0.0, 0.5};
static jj_state_t g_s = {0.0, 0.0, 0.0};
static double g_up[JJ_CURVE], g_down[JJ_CURVE];

QMC_EXPORT(qmc_jj_curve_len)
int qmc_jj_curve_len(void) { return JJ_CURVE; }

// Mean voltage on up sweep (i = imax k / (curve_len - 1)) and down sweep at
// same currents
QMC_EXPORT(qmc_jj_up_buffer)
double *qmc_jj_up_buffer(void) { return g_up; }
QMC_EXPORT(qmc_jj_down_buffer)
double *qmc_jj_down_buffer(void) { return g_down; }

// Set junction: quality factor in [0.05, 20], dc current in [-5, 5], ac
// amplitude in [0, 5], frequency in [0, 5]. 0 or -1
QMC_EXPORT(qmc_jj_set)
int qmc_jj_set(double q, double i_dc, double i_ac, double omega) {
  if (!(q >= 0.05 && q <= 20.0 && i_dc >= -5.0 && i_dc <= 5.0 && i_ac >= 0.0 &&
        i_ac <= 5.0 && omega >= 0.0 && omega <= 5.0)) {
    return -1;
  }

  g_p = (jj_params_t){q, i_dc, i_ac, omega};

  return 0;
}

QMC_EXPORT(qmc_jj_reset)
int qmc_jj_reset(double phi, double v) {
  if (!isfinite(phi) || !isfinite(v)) {
    return -1;
  }

  g_s = (jj_state_t){phi, v, 0.0};

  return 0;
}

// Advance live state by `time` (0 .. 10)
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_jj_advance)
int qmc_jj_advance(double time) {
  if (!(time >= 0.0 && time <= 10.0)) {
    return -1;
  }

  return jj_advance(&g_s, &g_p, time, 0.01);
}

QMC_EXPORT(qmc_jj_phi)
double qmc_jj_phi(void) { return g_s.phi; }
QMC_EXPORT(qmc_jj_v)
double qmc_jj_v(void) { return g_s.v; }
QMC_EXPORT(qmc_jj_time)
double qmc_jj_time(void) { return g_s.t; }

// I-V curves for currents 0 .. imax with current q, i_ac and \omega
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_jj_iv)
int qmc_jj_iv(double imax) {
  double tmp[JJ_CURVE];

  if (!(imax > 0.0 && imax <= 5.0)) {
    return -1;
  }

  jj_params_t p = g_p;
  p.i_dc = 0.0;

  if (jj_sweep(&p, 0.0, imax, JJ_CURVE, g_up) != 0 ||
      jj_sweep(&p, imax, 0.0, JJ_CURVE, tmp) != 0) {
    return -1;
  }

  for (int k = 0; k < JJ_CURVE; k++) {
    g_down[k] = tmp[JJ_CURVE - 1 - k]; // store by increasing current
  }

  return 0;
}
