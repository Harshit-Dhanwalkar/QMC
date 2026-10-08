/*
 * qmc_lz.c - WebAssembly front end for Landau-Zener / Stueckelberg demo
 *
 * Physics: physics/landau_zener.c (RK4 Bloch-vector integration under a swept
 * detuning, closed-form jump probability). This file keeps one running sweep
 * for animation and two scans (jump probability against sweep rate and
 * against sweep amplitude) that fill a scratch buffer
 */

#include "../physics/landau_zener.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { LZ_MAX_SCAN = 256 };

#define LZ_DT 0.01

static double g_v[3];
static double g_omega;
static double g_rate;
static double g_amp;
static int g_passes = 1;
static double g_t;
static double g_scan[LZ_MAX_SCAN];

static double pass_duration(void) { return 2.0 * g_amp / g_rate; }

QMC_EXPORT(qmc_lz_vec)
double *qmc_lz_vec(void) { return g_v; }

QMC_EXPORT(qmc_lz_scan_buffer)
double *qmc_lz_scan_buffer(void) { return g_scan; }

QMC_EXPORT(qmc_lz_scan_max)
int qmc_lz_scan_max(void) { return LZ_MAX_SCAN; }

QMC_EXPORT(qmc_lz_time)
double qmc_lz_time(void) { return g_t; }

QMC_EXPORT(qmc_lz_total_time)
double qmc_lz_total_time(void) { return g_passes * pass_duration(); }

// Detuning at the current time (zigzag between -amp and +amp)
QMC_EXPORT(qmc_lz_delta)
double qmc_lz_delta(void) {
  double dur = pass_duration();
  int p = (int)floor(g_t / dur);

  if (p >= g_passes) {
    p = g_passes - 1;
  }

  double local = g_t - p * dur;
  double dir = (p % 2 == 0) ? 1.0 : -1.0;

  return -dir * g_amp + dir * g_rate * local;
}

// Population of upper adiabatic level at current time
QMC_EXPORT(qmc_lz_upper)
double qmc_lz_upper(void) {
  return lz_upper_population(g_v, g_omega, qmc_lz_delta());
}

// Landau-Zener probability for current coupling and rate
QMC_EXPORT(qmc_lz_exact)
double qmc_lz_exact(void) { return lz_probability(g_omega, g_rate); }

/* Start a sweep in lower adiabatic state at delta = -amp. \omega in
 * [0, 3], rate in (0, 20], amp in [2, 20], passes 1 or 2. Returns 0 */
QMC_EXPORT(qmc_lz_reset)
int qmc_lz_reset(double omega, double rate, double amp, int passes) {
  if (!(omega >= 0.0 && omega <= 3.0) || !(rate > 0.0 && rate <= 20.0) ||
      !(amp >= 2.0 && amp <= 20.0) || passes < 1 || passes > 2) {
    return -1;
  }

  g_omega = omega;
  g_rate = rate;
  g_amp = amp;
  g_passes = passes;
  g_t = 0.0;

  return lz_lower_state(g_v, omega, -amp);
}

/* Advance by `duration` of sweep time (clipped to end)
 * Returns 0 while running, 1 once sweep has finished, -1 on error */
QMC_EXPORT(qmc_lz_step)
int qmc_lz_step(double duration) {
  double total = qmc_lz_total_time();
  double dur = pass_duration();

  if (!(duration > 0.0) || !(total > 0.0)) {
    return -1;
  }

  double left = fmin(duration, total - g_t);
  while (left > 1e-12) {
    int p = (int)floor((g_t + 1e-12) / dur);

    if (p >= g_passes) {
      break;
    }

    double pass_end = (p + 1) * dur;
    double chunk = fmin(left, pass_end - g_t);
    double dir = (p % 2 == 0) ? 1.0 : -1.0;
    double d0 = qmc_lz_delta();

    if (lz_advance(g_v, g_omega, d0, dir * g_rate, chunk, LZ_DT) != 0) {
      return -1;
    }

    g_t += chunk;
    left -= chunk;
  }

  return g_t >= total - 1e-9 ? 1 : 0;
}

// P(jump) after one pass for n log-spaced rates in [rmin, rmax]
QMC_EXPORT(qmc_lz_scan_rate)
int qmc_lz_scan_rate(double omega, double amp, int n, double rmin,
                     double rmax) {
  if (n < 2 || n > LZ_MAX_SCAN || !(rmin > 0.0) || !(rmax > rmin) ||
      !(omega >= 0.0 && omega <= 3.0) || !(amp >= 2.0 && amp <= 20.0)) {
    return -1;
  }

  for (int i = 0; i < n; i++) {
    double rate = rmin * pow(rmax / rmin, (double)i / (n - 1));
    double v[3];

    lz_lower_state(v, omega, -amp);
    if (lz_sweep(v, omega, rate, amp, LZ_DT, 1) != 0) {
      return -1;
    }

    g_scan[i] = lz_upper_population(v, omega, amp);
  }

  return 0;
}

// P(upper level) after two passes for n amplitudes in [a0, a1]
QMC_EXPORT(qmc_lz_scan_amp)
int qmc_lz_scan_amp(double omega, double rate, int n, double a0, double a1) {
  if (n < 2 || n > LZ_MAX_SCAN || !(a0 >= 2.0) || !(a1 > a0) || a1 > 20.0 ||
      !(omega >= 0.0 && omega <= 3.0) || !(rate > 0.0 && rate <= 20.0)) {
    return -1;
  }

  for (int i = 0; i < n; i++) {
    double amp = a0 + (a1 - a0) * i / (n - 1);
    double v[3];

    lz_lower_state(v, omega, -amp);
    if (lz_sweep(v, omega, rate, amp, LZ_DT, 2) != 0) {
      return -1;
    }

    g_scan[i] = lz_upper_population(v, omega, -amp);
  }

  return 0;
}
