/*
 * qmc_zeno.c - WebAssembly front end for quantum-Zeno demo
 *
 * Physics: physics/zeno.c (pulsed projective measurements, Bloch equations with
 * measurement dephasing, Monte-Carlo trajectories). Drive is fixed at \Omega =
 * 1 and total time at $T = \pi$, time of a complete flip, so every measurement
 * strategy is compared on same experiment
 */

#include "../physics/zeno.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { ZN_CURVE = 401, ZN_MAX_N = 400 };

static const double ZN_OMEGA = 1.0;
static const double ZN_T = 3.141592653589793;
static double g_time[ZN_CURVE], g_n_curve[ZN_CURVE], g_rate[ZN_CURVE],
    g_cont[ZN_CURVE];
static int g_alive[ZN_MAX_N + 1];
static double g_alive_f[ZN_MAX_N + 1];

QMC_EXPORT(qmc_zn_curve_len)
int qmc_zn_curve_len(void) { return ZN_CURVE; }

QMC_EXPORT(qmc_zn_max_n)
int qmc_zn_max_n(void) { return ZN_MAX_N; }

// Survival during run for N measurements: P(t), t = T i / (CURVE - 1)
QMC_EXPORT(qmc_zn_time_buffer)
double *qmc_zn_time_buffer(void) { return g_time; }

// Survival at t = T versus number of measurements, N = 1 .. CURVE - 1
// (index N - 1)
QMC_EXPORT(qmc_zn_n_buffer)
double *qmc_zn_n_buffer(void) { return g_n_curve; }

// Slowest decay rate (units of Omega) versus dephasing Gamma = gmax i / (CURVE
// - 1)
QMC_EXPORT(qmc_zn_rate_buffer)
double *qmc_zn_rate_buffer(void) { return g_rate; }

// Survival during run under continuous dephasing: P(t), t = T i / (CURVE - 1)
QMC_EXPORT(qmc_zn_cont_buffer)
double *qmc_zn_cont_buffer(void) { return g_cont; }

// Fraction of Monte-Carlo trajectories still alive after k measurements, k = 0
// .. N
QMC_EXPORT(qmc_zn_alive_buffer)
double *qmc_zn_alive_buffer(void) { return g_alive_f; }

QMC_EXPORT(qmc_zn_total_time)
double qmc_zn_total_time(void) { return ZN_T; }

// Pulsed survival at t = T with N measurements
QMC_EXPORT(qmc_zn_pulsed)
double qmc_zn_pulsed(int n) { return zeno_survival_pulsed(ZN_OMEGA, ZN_T, n); }

// Fill time curve for N measurements and all-N curve. Returns 0, or -1 on bad
// input
QMC_EXPORT(qmc_zn_pulsed_curves)
int qmc_zn_pulsed_curves(int n) {
  if (n < 1 || n > ZN_MAX_N) {
    return -1;
  }

  for (int i = 0; i < ZN_CURVE; i++) {
    g_time[i] =
        zeno_survival_pulsed_at(ZN_OMEGA, ZN_T, n, ZN_T * i / (ZN_CURVE - 1));
  }

  for (int k = 1; k < ZN_CURVE; k++) {
    g_n_curve[k - 1] = zeno_survival_pulsed(ZN_OMEGA, ZN_T, k);
  }

  return 0;
}

// Continuous measurement of strength \gamma (units of \Omega): fills time
// curve. Returns survival at T
QMC_EXPORT(qmc_zn_continuous)
double qmc_zn_continuous(double gamma) {
  if (!(gamma >= 0.0 && gamma <= 100.0)) {
    return NAN;
  }

  for (int i = 0; i < ZN_CURVE; i++) {
    g_cont[i] =
        zeno_survival_continuous(ZN_OMEGA, gamma, ZN_T * i / (ZN_CURVE - 1));
  }

  return g_cont[ZN_CURVE - 1];
}

// Rate curve over \Gamma in [0, gmax]
QMC_EXPORT(qmc_zn_rates)
int qmc_zn_rates(double gmax) {
  if (!(gmax > 0.0 && gmax <= 1000.0)) {
    return -1;
  }

  for (int i = 0; i < ZN_CURVE; i++) {
    g_rate[i] = zeno_continuous_rate(ZN_OMEGA, gmax * i / (ZN_CURVE - 1));
  }

  return 0;
}

QMC_EXPORT(qmc_zn_rate)
double qmc_zn_rate(double gamma) {
  return zeno_continuous_rate(ZN_OMEGA, gamma);
}

// Simulate n_traj qubits through N measurements
// Returns final alive fraction, NaN on bad input
QMC_EXPORT(qmc_zn_sample)
double qmc_zn_sample(int n, int n_traj, int seed) {
  if (n < 1 || n > ZN_MAX_N ||
      zeno_sample_pulsed(ZN_OMEGA, ZN_T, n, n_traj, (unsigned)seed, g_alive) !=
          0) {
    return NAN;
  }

  for (int k = 0; k <= n; k++) {
    g_alive_f[k] = (double)g_alive[k] / n_traj;
  }

  return g_alive_f[n];
}
