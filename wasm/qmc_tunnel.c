/*
 * qmc_tunnel.c - WebAssembly front end for resonant-tunnelling demo
 *
 * Physics: physics/double_barrier.c (transfer-matrix transmission through
 * square barriers). Structure is barrier 1 (height v1, thickness b1), a
 * well of width w at potential vw and barrier 2 (v2, b2), with \hbar^2 / 2m = 1
 * so E = k^2
 */

#include "../physics/double_barrier.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { DT_CURVE = 401 };

static double g_v[3] = {30.0, 0.0, 30.0};
static double g_d[3] = {0.5, 3.0, 0.5};
static double g_te[DT_CURVE];
static double g_dwell[DT_CURVE];
static double g_psi[DT_CURVE];
static double g_pot[DT_CURVE], g_r;

QMC_EXPORT(qmc_dt_curve_len)
int qmc_dt_curve_len(void) { return DT_CURVE; }

// Transmission at E = emax (i + 1) / curve_len and dwell time there
QMC_EXPORT(qmc_dt_te_buffer)
double *qmc_dt_te_buffer(void) { return g_te; }
QMC_EXPORT(qmc_dt_dwell_buffer)
double *qmc_dt_dwell_buffer(void) { return g_dwell; }
// |\psi|^2 and potential on x = xmin .. xmax (qmc_dt_profile)
QMC_EXPORT(qmc_dt_psi_buffer)
double *qmc_dt_psi_buffer(void) { return g_psi; }
QMC_EXPORT(qmc_dt_pot_buffer)
double *qmc_dt_pot_buffer(void) { return g_pot; }

// Define structure
// Returns 0, or -1 on bad input (heights >= 0 thicknesses in (0, 20], well
// width in (0, 20])
QMC_EXPORT(qmc_dt_set)
int qmc_dt_set(double v1, double b1, double w, double b2, double v2,
               double vw) {
  if (!(v1 >= 0.0 && v1 <= 500.0 && v2 >= 0.0 && v2 <= 500.0 && vw >= -50.0 &&
        vw <= 500.0 && b1 > 0.0 && b1 <= 20.0 && b2 > 0.0 && b2 <= 20.0 &&
        w > 0.0 && w <= 20.0)) {
    return -1;
  }

  g_v[0] = v1;
  g_v[1] = vw;
  g_v[2] = v2;
  g_d[0] = b1;
  g_d[1] = w;
  g_d[2] = b2;

  return 0;
}

QMC_EXPORT(qmc_dt_length)
double qmc_dt_length(void) { return g_d[0] + g_d[1] + g_d[2]; }

// n-th level of infinite well of current width (above its floor)
QMC_EXPORT(qmc_dt_well_level)
double qmc_dt_well_level(int n) {
  return n >= 1 ? g_v[1] + pow(n * 3.141592653589793 / g_d[1], 2) : NAN;
}

// T(E) and dwell time on E in (0, emax]
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_dt_curves)
int qmc_dt_curves(double emax) {
  if (!(emax > 0.0 && emax <= 1000.0)) {
    return -1;
  }

  for (int i = 0; i < DT_CURVE; i++) {
    double e = emax * (i + 1) / DT_CURVE;

    g_te[i] = db_transmission(3, g_v, g_d, e);
    g_dwell[i] = db_dwell_time(3, g_v, g_d, e);
  }

  return 0;
}

// Transmission at one energy (reflection via qmc_dt_last_reflection)
QMC_EXPORT(qmc_dt_point)
double qmc_dt_point(double e) {
  db_result_t r;

  if (db_solve(3, g_v, g_d, e, &r) != 0) {
    return NAN;
  }

  g_r = r.reflection;

  return r.transmission;
}

QMC_EXPORT(qmc_dt_last_reflection)
double qmc_dt_last_reflection(void) { return g_r; }

QMC_EXPORT(qmc_dt_dwell)
double qmc_dt_dwell(double e) { return db_dwell_time(3, g_v, g_d, e); }

// |\psi|^2 for unit incident amplitude and potential on [xmin, xmax]
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_dt_profile)
int qmc_dt_profile(double e, double xmin, double xmax) {
  double len = qmc_dt_length();

  if (db_profile(3, g_v, g_d, e, xmin, xmax, DT_CURVE, g_psi) != 0) {
    return -1;
  }

  for (int i = 0; i < DT_CURVE; i++) {
    double x = xmin + (xmax - xmin) * i / (DT_CURVE - 1);

    g_pot[i] =
        (x < 0.0 || x > len)
            ? 0.0
            : (x < g_d[0] ? g_v[0] : (x < g_d[0] + g_d[1] ? g_v[1] : g_v[2]));
  }

  return 0;
}
