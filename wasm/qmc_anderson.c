/*
 * qmc_anderson.c - WebAssembly front end for Anderson-localisation demo
 *
 * Physics: physics/anderson.c (disorder, RK4 propagation of a tight-binding
 * chain, width / IPR / energy, weak-disorder Lyapunov exponent). This file
 * only owns one chain of AN_N sites (hop t = 1) with packet started on middle
 * site, and scratch buffer JavaScript reads from
 */

#include "../physics/anderson.h"
#include "complex.h"
#include "vector.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { AN_N = 801 };

#define AN_DT 0.01
#define AN_HOP 1.0

static cvector_t *g_psi;
static double g_eps[AN_N];
static double g_buf[2 * AN_N]; /* |\psi|^2, then on-site energies */
static double g_t;

QMC_EXPORT(qmc_an_n)
int qmc_an_n(void) { return AN_N; }

QMC_EXPORT(qmc_an_dt)
double qmc_an_dt(void) { return AN_DT; }

QMC_EXPORT(qmc_an_time)
double qmc_an_time(void) { return g_t; }

/* Scratch buffer: [0,N) |psi|^2, [N,2N) on-site energies */
QMC_EXPORT(qmc_an_buffer)
double *qmc_an_buffer(void) { return g_buf; }

/* Refresh scratch buffer from current state */
QMC_EXPORT(qmc_an_density)
void qmc_an_density(void) {
  if (!g_psi) {
    return;
  }

  for (int j = 0; j < AN_N; j++) {
    g_buf[j] = c_abs2(g_psi->data[j]);
    g_buf[AN_N + j] = g_eps[j];
  }
}

/* New disorder (strength w in [0, 12], integer seed) and a packet on the
 * middle site at t = 0. Returns 0 on success */
QMC_EXPORT(qmc_an_reset)
int qmc_an_reset(double w, int seed) {
  if (!(w >= 0.0 && w <= 12.0)) {
    return -1;
  }
  if (!g_psi) {
    g_psi = cvector_alloc(AN_N);
  }
  if (!g_psi) {
    return -1;
  }
  if (anderson_disorder(g_eps, AN_N, w, (unsigned long long)(unsigned)seed) !=
      0) {
    return -1;
  }

  cvector_fill(g_psi, c_new(0.0, 0.0));
  g_psi->data[AN_N / 2] = c_new(1.0, 0.0);
  g_t = 0.0;
  qmc_an_density();

  return 0;
}

/* Advance n steps of AN_DT. Returns 0 on success */
QMC_EXPORT(qmc_an_step)
int qmc_an_step(int n) {
  if (!g_psi) {
    return -1;
  }

  int rc = anderson_evolve(g_psi, g_eps, AN_HOP, AN_DT, n);
  if (rc == 0) {
    g_t += AN_DT * n;
    qmc_an_density();
  }

  return rc;
}

QMC_EXPORT(qmc_an_norm)
double qmc_an_norm(void) { return anderson_norm(g_psi); }

QMC_EXPORT(qmc_an_width)
double qmc_an_width(void) { return anderson_width(g_psi); }

QMC_EXPORT(qmc_an_ipr)
double qmc_an_ipr(void) { return anderson_ipr(g_psi); }

QMC_EXPORT(qmc_an_energy)
double qmc_an_energy(void) { return anderson_energy(g_psi, g_eps, AN_HOP); }

/* Weak-disorder localisation length (1 / \gamma) at energy e (NaN outside band
 * or for w = 0, where it is infinite) */
QMC_EXPORT(qmc_an_xi)
double qmc_an_xi(double w, double e) {
  double gamma = anderson_lyapunov_weak(w, AN_HOP, e);

  return gamma > 0.0 ? 1.0 / gamma : NAN;
}
