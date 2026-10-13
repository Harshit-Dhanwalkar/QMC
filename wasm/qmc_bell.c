/*
 * qmc_bell.c - WebAssembly front end for CHSH Bell-test demo
 *
 * Physics: physics/bell_chsh.c (Werner-state correlations from density matrix,
 * CHSH combination, a local hidden-angle model and a pair sampler)
 * This file owns 4 measurement angles, visibility, correlation curves and last
 * simulated sample
 */

#include "../physics/bell_chsh.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { BL_CURVE_MAX = 361 };

static double g_ang[4] = {0.0, 1.5707963267948966, 0.7853981633974483,
                          2.356194490192345};
static double g_v = 1.0;
static double g_q[BL_CURVE_MAX], g_l[BL_CURVE_MAX];
static double g_est[4], g_sigma, g_s_est;

QMC_EXPORT(qmc_bl_curve_max)
int qmc_bl_curve_max(void) { return BL_CURVE_MAX; }

// Correlation versus angle difference d in [0, 2 \pi]: quantum (visibility V)
// and hidden-angle local model
QMC_EXPORT(qmc_bl_quantum_buffer)
double *qmc_bl_quantum_buffer(void) { return g_q; }

QMC_EXPORT(qmc_bl_local_buffer)
double *qmc_bl_local_buffer(void) { return g_l; }

// Angles are in radians: Alice a, a2; Bob b, b2. Returns 0 on success, -1 on
// bad input
QMC_EXPORT(qmc_bl_set)
int qmc_bl_set(double a, double a2, double b, double b2, double visibility) {
  if (!isfinite(a) || !isfinite(a2) || !isfinite(b) || !isfinite(b2) ||
      !(visibility >= 0.0 && visibility <= 1.0)) {
    return -1;
  }

  g_ang[0] = a;
  g_ang[1] = a2;
  g_ang[2] = b;
  g_ang[3] = b2;
  g_v = visibility;

  return 0;
}

// S with current settings (quantum / local model)
QMC_EXPORT(qmc_bl_chsh)
double qmc_bl_chsh(void) {
  return bell_chsh(g_v, g_ang[0], g_ang[1], g_ang[2], g_ang[3]);
}

QMC_EXPORT(qmc_bl_local_chsh)
double qmc_bl_local_chsh(void) {
  return bell_local_chsh(g_ang[0], g_ang[1], g_ang[2], g_ang[3]);
}

// Maximum |S| over all settings on a \pi/4 grid for current visibility:
// Tsirelson, 2 \sqrt{2} V
QMC_EXPORT(qmc_bl_max)
double qmc_bl_max(void) { return bell_max_chsh(g_v, 8); }

// Correlation for setting pair i = 0..3: (a,b), (a,b2), (a2,b), (a2,b2)
QMC_EXPORT(qmc_bl_e)
double qmc_bl_e(int i) {
  static const int pa[4] = {0, 0, 1, 1}, pb[4] = {2, 3, 2, 3};

  if (i < 0 || i > 3) {
    return NAN;
  }

  return bell_correlation(g_v, g_ang[pa[i]], g_ang[pb[i]]);
}

QMC_EXPORT(qmc_bl_e_local)
double qmc_bl_e_local(int i) {
  static const int pa[4] = {0, 0, 1, 1}, pb[4] = {2, 3, 2, 3};

  if (i < 0 || i > 3) {
    return NAN;
  }

  return bell_local_correlation(g_ang[pa[i]], g_ang[pb[i]]);
}

QMC_EXPORT(qmc_bl_curve)
int qmc_bl_curve(int count) {
  if (count < 2 || count > BL_CURVE_MAX) {
    return -1;
  }

  for (int i = 0; i < count; i++) {
    double d = 6.283185307179586 * i / (count - 1);

    g_q[i] = bell_correlation(g_v, 0.0, d);
    g_l[i] = bell_local_correlation(0.0, d);
  }

  return 0;
}

// Simulate n pairs per setting; returns estimate of S (sign included)
QMC_EXPORT(qmc_bl_sample)
double qmc_bl_sample(int n, int seed) {
  double sums[4];
  double var = 0.0;

  g_s_est = bell_sample(g_v, g_ang[0], g_ang[1], g_ang[2], g_ang[3], n,
                        (unsigned)seed, sums);
  if (isnan(g_s_est)) {
    return NAN;
  }

  for (int i = 0; i < 4; i++) {
    g_est[i] = sums[i] / n;
    var += (1.0 - qmc_bl_e(i) * qmc_bl_e(i)) / n;
  }

  g_sigma = sqrt(var);

  return g_s_est;
}

QMC_EXPORT(qmc_bl_est)
double qmc_bl_est(int i) { return i < 0 || i > 3 ? NAN : g_est[i]; }

// One-sigma statistical error of sampled S
QMC_EXPORT(qmc_bl_sigma)
double qmc_bl_sigma(void) { return g_sigma; }
