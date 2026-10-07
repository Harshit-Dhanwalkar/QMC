/*
 * qmc_qwalk.c - WebAssembly front end for the quantum-walk demo
 *
 * Physics: physics/quantum_walk.c (coin + shift step, observables, classical
 * binomial reference). This file only owns one walk of QW_N sites started at
 * the centre and the scratch buffer JavaScript reads from
 */

#include "../physics/quantum_walk.h"
#include "vector.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { QW_N = 1024, QW_TMAX = 480 };

static cvector_t *g_up;
static cvector_t *g_down;
static double g_buf[2 * QW_N]; /* quantum P(x), classical P(x) */
static double g_theta = 0.7853981633974483;
static int g_t;

QMC_EXPORT(qmc_qw_n)
int qmc_qw_n(void) { return QW_N; }

QMC_EXPORT(qmc_qw_tmax)
int qmc_qw_tmax(void) { return QW_TMAX; }

QMC_EXPORT(qmc_qw_time)
int qmc_qw_time(void) { return g_t; }

/* Scratch buffer: [0,N) quantum P(x), [N,2N) classical P(x) */
QMC_EXPORT(qmc_qw_buffer)
double *qmc_qw_buffer(void) { return g_buf; }

static int refresh(void) {
  return qwalk_probability(g_up, g_down, g_buf) +
         qwalk_classical(g_buf + QW_N, QW_N, g_t);
}

/* Restart with coin angle theta in [0, pi/2] and start state 0 (up), 1 (down)
 * or 2 (symmetric). Returns 0 on success */
QMC_EXPORT(qmc_qw_reset)
int qmc_qw_reset(double theta, int coin) {
  if (!(theta >= 0.0 && theta <= 1.5707963267948966)) {
    return -1;
  }
  if (!g_up) {
    g_up = cvector_alloc(QW_N);
  }
  if (!g_down) {
    g_down = cvector_alloc(QW_N);
  }
  if (!g_up || !g_down || qwalk_init(g_up, g_down, coin) != 0) {
    return -1;
  }

  g_theta = theta;
  g_t = 0;

  return refresh() == 0 ? 0 : -1;
}

/* Advance n steps (total time is capped at QW_TMAX). Returns 0 on success */
QMC_EXPORT(qmc_qw_step)
int qmc_qw_step(int n) {
  if (!g_up || n < 1 || g_t + n > QW_TMAX) {
    return -1;
  }
  if (qwalk_step(g_up, g_down, g_theta, n) != 0) {
    return -1;
  }
  g_t += n;

  return refresh() == 0 ? 0 : -1;
}

QMC_EXPORT(qmc_qw_norm)
double qmc_qw_norm(void) { return qwalk_norm(g_up, g_down); }

QMC_EXPORT(qmc_qw_mean)
double qmc_qw_mean(void) { return qwalk_mean(g_up, g_down); }

QMC_EXPORT(qmc_qw_variance)
double qmc_qw_variance(void) { return qwalk_variance(g_up, g_down); }

/* Limit of variance / t^2 for the current coin angle */
QMC_EXPORT(qmc_qw_asymptote)
double qmc_qw_asymptote(void) { return qwalk_asymptotic_variance(g_theta); }
