/*
 * qmc_tfim.c - WebAssembly front end for the Ising-quench demo.
 *
 * Physics: physics/tfim_quench.c (free-fermion solution of the sudden quench of
 * the transverse-field Ising chain). This file owns the two field values, a
 * 256-site ring for correlations and magnetisation, a 1024-site ring for the
 * smooth Loschmidt-rate curve, and the buffers JavaScript reads.
 */

#include "../physics/tfim_quench.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { TQ_N = 256, TQ_CURVE_N = 1024, TQ_R = 40, TQ_CURVE_MAX = 400 };

static double g_hi = 0.3;
static double g_hf = 2.0;
static double g_zz[TQ_R + 1];
static double g_rate_curve[TQ_CURVE_MAX];
static double g_mx_curve[TQ_CURVE_MAX];
static double g_density[TQ_N / 2];
static double g_mx;
static double g_rate;

QMC_EXPORT(qmc_tq_rmax)
int qmc_tq_rmax(void) { return TQ_R; }

QMC_EXPORT(qmc_tq_curve_max)
int qmc_tq_curve_max(void) { return TQ_CURVE_MAX; }

QMC_EXPORT(qmc_tq_modes)
int qmc_tq_modes(void) { return TQ_N / 2; }

/* <sz_0 sz_r> for r = 0 .. rmax at the last evaluated time. */
QMC_EXPORT(qmc_tq_zz_buffer)
double *qmc_tq_zz_buffer(void) { return g_zz; }

QMC_EXPORT(qmc_tq_rate_buffer)
double *qmc_tq_rate_buffer(void) { return g_rate_curve; }

QMC_EXPORT(qmc_tq_mx_buffer)
double *qmc_tq_mx_buffer(void) { return g_mx_curve; }

/* Quasiparticle density (1 - n_i . n_f) / 2 of mode m, k = (2m + 1) pi / N. */
QMC_EXPORT(qmc_tq_density_buffer)
double *qmc_tq_density_buffer(void) { return g_density; }

/* Set the fields, J = 1: h_i and h_f in [0, 3]. Returns 0 on success. */
QMC_EXPORT(qmc_tq_set)
int qmc_tq_set(double h_i, double h_f) {
  if (!(h_i >= 0.0 && h_i <= 3.0) || !(h_f >= 0.0 && h_f <= 3.0)) {
    return -1;
  }
  g_hi = h_i;
  g_hf = h_f;
  for (int m = 0; m < TQ_N / 2; m++) {
    double k = (2 * m + 1) * 3.14159265358979323846 / TQ_N;
    double yi = sin(k), zi = h_i - cos(k), yf = sin(k), zf = h_f - cos(k);
    double c = (yi * yf + zi * zf) / (sqrt(yi * yi + zi * zi) * sqrt(yf * yf + zf * zf));

    g_density[m] = 0.5 * (1.0 - c);
  }

  return 0;
}

/* Evaluate the state at time t >= 0: fills the correlation buffer and
 * qmc_tq_last_mx / qmc_tq_last_rate. Returns 0 on success. */
QMC_EXPORT(qmc_tq_eval)
int qmc_tq_eval(double t) {
  if (!(t >= 0.0 && t <= 1000.0) ||
      tfim_quench_zz(1.0, g_hi, g_hf, t, TQ_N, TQ_R, g_zz) != 0) {
    return -1;
  }
  g_mx = tfim_quench_mx(1.0, g_hi, g_hf, t, TQ_N);
  g_rate = tfim_quench_loschmidt_rate(1.0, g_hi, g_hf, t, TQ_CURVE_N);

  return 0;
}

QMC_EXPORT(qmc_tq_last_mx)
double qmc_tq_last_mx(void) { return g_mx; }

QMC_EXPORT(qmc_tq_last_rate)
double qmc_tq_last_rate(void) { return g_rate; }

/* Rate function and <sx> at `count` times from 0 to t_max, into the curve
 * buffers. Returns 0 on success. */
QMC_EXPORT(qmc_tq_curve)
int qmc_tq_curve(double t_max, int count) {
  if (count < 2 || count > TQ_CURVE_MAX || !(t_max > 0.0 && t_max <= 1000.0)) {
    return -1;
  }
  for (int i = 0; i < count; i++) {
    double t = t_max * i / (count - 1);

    g_rate_curve[i] = tfim_quench_loschmidt_rate(1.0, g_hi, g_hf, t, TQ_CURVE_N);
    g_mx_curve[i] = tfim_quench_mx(1.0, g_hi, g_hf, t, TQ_CURVE_N);
  }

  return 0;
}

/* Critical time t*_n of the dynamical transition (NaN without a crossing). */
QMC_EXPORT(qmc_tq_tcrit)
double qmc_tq_tcrit(int n) { return tfim_quench_critical_time(1.0, g_hi, g_hf, n); }

/* Front speed of the correlation light cone, 4 min(J, h_f). */
QMC_EXPORT(qmc_tq_front_speed)
double qmc_tq_front_speed(void) { return 4.0 * (g_hf < 1.0 ? g_hf : 1.0); }
