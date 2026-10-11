/*
 * qmc_haldane.c - WebAssembly front end for Haldane Chern-insulator demo
 *
 * Physics: physics/tight_binding.c (Haldane model, ribbon spectra, Chern
 * number). This file owns one Haldane model (t1 = 1) and a ribbon of HD_W unit
 * cells, plus buffers JavaScript reads
 */

#include "../physics/tight_binding.h"
#include <math.h>
#include <stddef.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { HD_W = 14, HD_S = 2 * HD_W, HD_NK_MAX = 121 };

static tb_model_t *g_model;
static double g_t2 = 0.1;
static double g_phi = M_PI / 2;
static double g_mass;
static double g_vmax;
static double g_cmin;
static double g_energy[HD_NK_MAX * HD_S];
static double g_mean[HD_NK_MAX * HD_S];
static double g_prob[HD_W];

QMC_EXPORT(qmc_hd_width)
int qmc_hd_width(void) { return HD_W; }

QMC_EXPORT(qmc_hd_states)
int qmc_hd_states(void) { return HD_S; }

QMC_EXPORT(qmc_hd_nk_max)
int qmc_hd_nk_max(void) { return HD_NK_MAX; }

// Ribbon energies, state s at wave number index ik: [ik * states + s]
QMC_EXPORT(qmc_hd_energy_buffer)
double *qmc_hd_energy_buffer(void) { return g_energy; }

// Mean cell of each state in [-1, 1] (-1 first edge, +1 far edge)
QMC_EXPORT(qmc_hd_mean_buffer)
double *qmc_hd_mean_buffer(void) { return g_mean; }

// Probability per cell of state chosen by qmc_hd_profile
QMC_EXPORT(qmc_hd_prob_buffer)
double *qmc_hd_prob_buffer(void) { return g_prob; }

/* Set t2 in [0, 1/3], \phi in [-\pi, \pi] and sublattice mass M in [-3, 3]
 * (t1 = 1). Returns 0 on success */
QMC_EXPORT(qmc_hd_set)
int qmc_hd_set(double t2, double phi, double mass) {
  if (!(t2 >= 0.0 && t2 <= 1.0 / 3.0 + 1e-12) || !(fabs(phi) <= M_PI + 1e-12) ||
      !(fabs(mass) <= 3.0)) {
    return -1;
  }

  tb_model_t *m = tb_model_haldane(1.0, t2, phi, mass);
  if (!m) {
    return -1;
  }

  tb_model_free(g_model);

  g_model = m;
  g_t2 = t2;
  g_phi = phi;
  g_mass = mass;
  g_vmax = -1e30;
  g_cmin = 1e30;

  for (int i = 0; i < 24; i++) {
    for (int j = 0; j < 24; j++) {
      double e[2];

      tb_bands(g_model, 2.0 * M_PI * i / 24, 2.0 * M_PI * j / 24, e);

      g_vmax = fmax(g_vmax, e[0]);
      g_cmin = fmin(g_cmin, e[1]);
    }
  }

  return 0;
}

// Ribbon bands at n_k wave numbers k1 = 2 pi i / (n_k - 1). Returns 0
QMC_EXPORT(qmc_hd_compute)
int qmc_hd_compute(int n_k) {
  if (!g_model || n_k < 2 || n_k > HD_NK_MAX) {
    return -1;
  }

  for (int i = 0; i < n_k; i++) {
    if (tb_ribbon_bands(g_model, 2.0 * M_PI * i / (n_k - 1), HD_W,
                        g_energy + (size_t)i * HD_S,
                        g_mean + (size_t)i * HD_S) != TB_OK) {
      return -1;
    }
  }

  return 0;
}

// Cell-resolved probability of state `state` at wave-number index ik of n_k
QMC_EXPORT(qmc_hd_profile)
int qmc_hd_profile(int ik, int n_k, int state) {
  if (!g_model || n_k < 2 || ik < 0 || ik >= n_k) {
    return -1;
  }

  return tb_ribbon_profile(g_model, 2.0 * M_PI * ik / (n_k - 1), HD_W, state,
                           g_prob) == TB_OK
             ? 0
             : -1;
}

// Chern number of lower band (INT_MIN on phase boundary)
QMC_EXPORT(qmc_hd_chern)
int qmc_hd_chern(void) { return tb_haldane_chern(1.0, g_t2, g_phi, g_mass); }

// Phase boundary |M| = 3 \sqrt(3) t2 |\sin(\phi)|
QMC_EXPORT(qmc_hd_limit)
double qmc_hd_limit(void) { return 3.0 * sqrt(3.0) * g_t2 * fabs(sin(g_phi)); }

// Gap at valley K (+1) or K' (-1)
QMC_EXPORT(qmc_hd_valley_gap)
double qmc_hd_valley_gap(int valley) {
  return tb_haldane_valley_gap(g_t2, g_phi, g_mass, valley);
}

// Bulk valence-band maximum and conduction-band minimum
QMC_EXPORT(qmc_hd_vmax)
double qmc_hd_vmax(void) { return g_vmax; }

QMC_EXPORT(qmc_hd_cmin)
double qmc_hd_cmin(void) { return g_cmin; }
