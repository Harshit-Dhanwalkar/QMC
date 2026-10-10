/*
 * qmc_rotor.c - WebAssembly front end for kicked-rotor demo
 *
 * Physics: physics/kicked_rotor.c (FFT split-step quantum map and classical
 * standard map). This file owns one quantum state of KR_N angle points and one
 * classical ensemble of KR_NP points
 */

#include "../physics/kicked_rotor.h"
#include "vector.h"
#include <math.h>
#include <stddef.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { KR_N = 2048, KR_NP = 3000, KR_TMAX = 400 };

static cvector_t *g_psi;
static double g_theta[KR_NP];
static double g_p[KR_NP];
static double g_qprob[KR_N];
static double g_cbuf[2 * KR_NP]; /* \theta, p pairs */
static double g_k;
static double g_hbar = 1.0;
static int g_t;

QMC_EXPORT(qmc_kr_n)
int qmc_kr_n(void) { return KR_N; }

QMC_EXPORT(qmc_kr_np)
int qmc_kr_np(void) { return KR_NP; }

QMC_EXPORT(qmc_kr_tmax)
int qmc_kr_tmax(void) { return KR_TMAX; }

QMC_EXPORT(qmc_kr_time)
int qmc_kr_time(void) { return g_t; }

/* Quantum momentum probabilities, index j is m = j - N/2 */
QMC_EXPORT(qmc_kr_qbuffer)
double *qmc_kr_qbuffer(void) { return g_qprob; }

/* Classical ensemble as (\theta, p) pairs; p is not wrapped */
QMC_EXPORT(qmc_kr_cbuffer)
double *qmc_kr_cbuffer(void) { return g_cbuf; }

static int refresh(void) {
  for (int i = 0; i < KR_NP; i++) {
    g_cbuf[(ptrdiff_t)2 * i] = g_theta[i];
    g_cbuf[2 * i + 1] = g_p[i];
  }

  return krotor_momentum_distribution(g_psi, g_qprob);
}

/* Restart with kick strength k in [0, 12], effective Planck constant \hbar in
 * [0.5, 13] (4 \pi = 12.566 is quantum resonance) and a seed for classical
 * ensemble. Returns 0 on success */
QMC_EXPORT(qmc_kr_reset)
int qmc_kr_reset(double k, double hbar, int seed) {
  if (!(k >= 0.0 && k <= 12.0) || !(hbar >= 0.5 && hbar <= 13.0)) {
    return -1;
  }
  if (!g_psi) {
    g_psi = cvector_alloc(KR_N);
  }
  if (!g_psi || krotor_init(g_psi) != 0 ||
      krotor_classical_init(g_theta, g_p, KR_NP, hbar,
                            (unsigned long long)(unsigned)seed) != 0) {
    return -1;
  }

  g_k = k;
  g_hbar = hbar;
  g_t = 0;

  return refresh();
}

/* Advance n kicks (total capped at KR_TMAX). Returns 0 on success */
QMC_EXPORT(qmc_kr_step)
int qmc_kr_step(int n) {
  if (!g_psi || n < 1 || g_t + n > KR_TMAX) {
    return -1;
  }
  if (krotor_step(g_psi, g_k, g_hbar, n) != 0 ||
      krotor_classical_step(g_theta, g_p, KR_NP, g_k, n) != 0) {
    return -1;
  }

  g_t += n;

  return refresh();
}

/* <m^2> of quantum state */
QMC_EXPORT(qmc_kr_qm2)
double qmc_kr_qm2(void) { return krotor_momentum2(g_psi); }

/* Classical <p^2> / hbar^2, in same units as <m^2> */
QMC_EXPORT(qmc_kr_cm2)
double qmc_kr_cm2(void) {
  return krotor_classical_p2(g_p, KR_NP) / (g_hbar * g_hbar);
}

QMC_EXPORT(qmc_kr_norm)
double qmc_kr_norm(void) { return krotor_norm(g_psi); }

/* Quasilinear diffusion rate K^2 / (2 \hbar^2) of <m^2> per kick */
QMC_EXPORT(qmc_kr_diffusion)
double qmc_kr_diffusion(void) { return 0.5 * g_k * g_k / (g_hbar * g_hbar); }

/* Probability within 6% of edge of momentum grid (aliasing warning) */
QMC_EXPORT(qmc_kr_edge)
double qmc_kr_edge(void) {
  double sum = 0.0;
  int w = KR_N * 6 / 100;

  for (int j = 0; j < w; j++) {
    sum += g_qprob[j] + g_qprob[KR_N - 1 - j];
  }

  return sum;
}
