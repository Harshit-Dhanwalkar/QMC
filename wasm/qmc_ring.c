/*
 * qmc_ring.c - WebAssembly front end for  Aharonov-Bohm ring demo
 *
 * tight-binding ring with flux: spectrum, persistent current, two-terminal
 * transmission through leads at sites 0 and b)
 * Hopping is the unit of energy (t = 1)
 * Rings have 4 .. 32 sites
 */

#include "../physics/ab_ring.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { AR_CURVE = 201, AR_MAP = 64, AR_NMAX = 32 };

static double g_spec[AR_NMAX * AR_CURVE], g_cur[AR_CURVE], g_tf[AR_CURVE],
    g_te[AR_CURVE], g_map[AR_MAP * AR_MAP];
static double g_levels[AR_NMAX], g_r;

QMC_EXPORT(qmc_ar_curve_len)
int qmc_ar_curve_len(void) { return AR_CURVE; }
QMC_EXPORT(qmc_ar_map_n)
int qmc_ar_map_n(void) { return AR_MAP; }
QMC_EXPORT(qmc_ar_max_sites)
int qmc_ar_max_sites(void) { return AR_NMAX; }

// Sorted levels, level l at index l * curve_len + i, flux \phi = -1 + 2 i
// (curve_len - 1)
QMC_EXPORT(qmc_ar_spec_buffer)
double *qmc_ar_spec_buffer(void) { return g_spec; }
// Persistent current of ne electrons at the same fluxes
QMC_EXPORT(qmc_ar_cur_buffer)
double *qmc_ar_cur_buffer(void) { return g_cur; }
// Transmission against flux \phi = i / (curve_len - 1) in [0, 1]
QMC_EXPORT(qmc_ar_tf_buffer)
double *qmc_ar_tf_buffer(void) { return g_tf; }
// Transmission against energy E = -1.995 + 3.99 i / (curve_len - 1)
QMC_EXPORT(qmc_ar_te_buffer)
double *qmc_ar_te_buffer(void) { return g_te; }
// Transmission map, row = flux (0 .. 1), column = energy (-1.98 .. 1.98): index
// row * map_n + column
QMC_EXPORT(qmc_ar_map_buffer)
double *qmc_ar_map_buffer(void) { return g_map; }
// Sorted levels at one flux
QMC_EXPORT(qmc_ar_levels_buffer)
double *qmc_ar_levels_buffer(void) { return g_levels; }

static int ok(int n, int b) {
  return n >= 4 && n <= AR_NMAX && b >= 1 && b < n;
}

// Spectrum curves for an n-site ring
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_ar_spectrum)
int qmc_ar_spectrum(int n) {
  double e[AR_NMAX];

  if (n < 4 || n > AR_NMAX) {
    return -1;
  }

  for (int i = 0; i < AR_CURVE; i++) {
    ar_levels_sorted(n, -1.0 + 2.0 * i / (AR_CURVE - 1), 1.0, e);
    for (int l = 0; l < n; l++) {
      g_spec[l * AR_CURVE + i] = e[l];
    }
  }

  return 0;
}

// Persistent-current curve for ne electrons
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_ar_current)
int qmc_ar_current(int n, int ne) {
  if (n < 4 || n > AR_NMAX || ne < 0 || ne > n) {
    return -1;
  }

  for (int i = 0; i < AR_CURVE; i++) {
    // nudge off degeneracy at integer flux
    g_cur[i] = ar_current(n, ne, -1.0 + 2.0 * i / (AR_CURVE - 1) + 1e-9, 1.0);
  }

  return 0;
}

// Sorted levels at one flux; returns energy of ne lowest, NaN on bad input
QMC_EXPORT(qmc_ar_levels)
double qmc_ar_levels(int n, int ne, double flux) {
  if (n < 4 || n > AR_NMAX || ne < 0 || ne > n ||
      ar_levels_sorted(n, flux, 1.0, g_levels) != 0) {
    return NAN;
  }

  return ar_energy(n, ne, flux, 1.0);
}

// Current of ne electrons at one flux (nudged off exact degeneracies)
QMC_EXPORT(qmc_ar_current_at)
double qmc_ar_current_at(int n, int ne, double flux) {
  return (n < 4 || n > AR_NMAX) ? NAN : ar_current(n, ne, flux + 1e-9, 1.0);
}

// Transmission at one flux and energy; reflection probability is left for
// qmc_ar_last_reflection
QMC_EXPORT(qmc_ar_point)
double qmc_ar_point(int n, int b, double flux, double energy) {
  ar_scatter_t s;

  if (!ok(n, b) || ar_scatter(n, 0, b, flux, energy, 1.0, &s) != 0) {
    return NAN;
  }

  g_r = s.reflection;

  return s.transmission;
}

QMC_EXPORT(qmc_ar_last_reflection)
double qmc_ar_last_reflection(void) { return g_r; }

// T against flux at energy e
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_ar_trans_flux)
int qmc_ar_trans_flux(int n, int b, double e) {
  ar_scatter_t s;

  if (!ok(n, b) || !(fabs(e) < 2.0)) {
    return -1;
  }

  for (int i = 0; i < AR_CURVE; i++) {
    g_tf[i] = ar_scatter(n, 0, b, (double)i / (AR_CURVE - 1), e, 1.0, &s) == 0
                  ? s.transmission
                  : NAN;
  }

  return 0;
}

// T against energy at flux \phi
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_ar_trans_energy)
int qmc_ar_trans_energy(int n, int b, double flux) {
  ar_scatter_t s;

  if (!ok(n, b) || !isfinite(flux)) {
    return -1;
  }

  for (int i = 0; i < AR_CURVE; i++) {
    g_te[i] = ar_scatter(n, 0, b, flux, -1.995 + 3.99 * i / (AR_CURVE - 1), 1.0,
                         &s) == 0
                  ? s.transmission
                  : NAN;
  }

  return 0;
}

// Transmission map over flux and energy
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_ar_map)
int qmc_ar_map(int n, int b) {
  ar_scatter_t s;

  if (!ok(n, b)) {
    return -1;
  }

  for (int r = 0; r < AR_MAP; r++) {
    for (int c = 0; c < AR_MAP; c++) {
      g_map[r * AR_MAP + c] =
          ar_scatter(n, 0, b, (double)r / (AR_MAP - 1),
                     -1.98 + 3.96 * c / (AR_MAP - 1), 1.0, &s) == 0
              ? s.transmission
              : NAN;
    }
  }

  return 0;
}
