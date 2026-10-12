/*
 * qmc_wigner.c - WebAssembly front end for phase-space demo
 *
 * Physics: physics/phase_space.c (Fock-basis states, Kerr evolution, Wigner
 * function). This file owns one state, its Wigner grid and marginals page
 * draws. Units: \hbar = m = \omega = 1, x = (a + a^\dagger) / \sqrt(2)
 */

#include "../physics/phase_space.h"
#include "complex.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { WG_N = 128, WG_NMAX = 100, WG_FOCK_MAX = 12 };
enum {
  WG_COHERENT = 0,
  WG_FOCK = 1,
  WG_CAT_EVEN = 2,
  WG_CAT_ODD = 3,
  WG_SQUEEZED = 4
};

static const double WG_DX = 0.125;

static complex_t g_c[WG_NMAX + 1];
static complex_t g_psi[WG_N];
static double g_x[WG_N];
static double g_w[WG_N * WG_N];
static double g_mx[WG_N];
static double g_mp[WG_N];
static double g_pn[WG_NMAX + 1];
static double g_nbar, g_neg, g_wmin, g_wmax;

QMC_EXPORT(qmc_wg_n)
int qmc_wg_n(void) { return WG_N; }

QMC_EXPORT(qmc_wg_nmax)
int qmc_wg_nmax(void) { return WG_NMAX; }

// Half-width of phase-space window, in both x and p
QMC_EXPORT(qmc_wg_extent)
double qmc_wg_extent(void) { return (WG_N - 1) / 2.0 * WG_DX; }

// W[l * N + j] at p_l (row), x_j (column)
QMC_EXPORT(qmc_wg_w_buffer)
double *qmc_wg_w_buffer(void) { return g_w; }

// Marginals: |\psi(x)|^2 and momentum distribution, on same grid
QMC_EXPORT(qmc_wg_x_buffer)
double *qmc_wg_x_buffer(void) { return g_mx; }

QMC_EXPORT(qmc_wg_p_buffer)
double *qmc_wg_p_buffer(void) { return g_mp; }

// Photon-number distribution |c_n|^2, n = 0 .. nmax
QMC_EXPORT(qmc_wg_pn_buffer)
double *qmc_wg_pn_buffer(void) { return g_pn; }

QMC_EXPORT(qmc_wg_nbar)
double qmc_wg_nbar(void) { return g_nbar; }

// Negative volume (integral |W| - 1) / 2
QMC_EXPORT(qmc_wg_negativity)
double qmc_wg_negativity(void) { return g_neg; }

QMC_EXPORT(qmc_wg_wmin)
double qmc_wg_wmin(void) { return g_wmin; }

QMC_EXPORT(qmc_wg_wmax)
double qmc_wg_wmax(void) { return g_wmax; }

/*
 * Build a state, apply Kerr phase for \chi_t and compute its Wigner function
 *   kind 0 coherent: amp = |\alpha| (0..3.8), param = phase of alpha
 *   kind 1 Fock:     amp = n (0..12)
 *   kind 2/3 cat:    amp = |\alpha| (0..3.8), param = phase of alpha (even /
 *                    odd)
 *   kind 4 squeezed: amp = r (0..1.2), param = squeezing angle
 *
 * Returns 0 on success, -1 on invalid input
 */
QMC_EXPORT(qmc_wg_state)
int qmc_wg_state(int kind, double amp, double param, double chi_t) {
  int status = -1;

  if (!(amp >= 0.0) || !(isfinite(param)) || !(isfinite(chi_t))) {
    return -1;
  }

  switch (kind) {
  case WG_COHERENT:
    if (amp <= 3.8) {
      status =
          ps_coherent_state(g_c, WG_NMAX, amp * cos(param), amp * sin(param));
    }

    break;
  case WG_FOCK:
    if (amp <= WG_FOCK_MAX + 0.5) {
      status = ps_fock_state(g_c, WG_NMAX, (int)(amp + 0.5));
    }

    break;
  case WG_CAT_EVEN:
  case WG_CAT_ODD:
    if (amp <= 3.8) {
      status = ps_cat_state(g_c, WG_NMAX, amp * cos(param), amp * sin(param),
                            kind == WG_CAT_EVEN ? 1 : -1);
    }

    break;
  case WG_SQUEEZED:
    if (amp <= 1.2) {
      status = ps_squeezed_vacuum(g_c, WG_NMAX, amp, param);
    }

    break;

  default:
    break;
  }

  if (status != 0 || ps_kerr_evolve(g_c, WG_NMAX, chi_t) != 0) {
    return -1;
  }

  for (int j = 0; j < WG_N; j++) {
    g_x[j] = (j - (WG_N - 1) / 2.0) * WG_DX;
  }

  if (ps_wavefunction(g_c, WG_NMAX, g_x, WG_N, g_psi) != 0 ||
      ps_wigner(g_psi, WG_N, WG_DX, g_x, WG_N, g_w) != 0) {
    return -1;
  }

  g_wmin = g_wmax = g_w[0];
  for (int j = 0; j < WG_N; j++) {
    g_mx[j] = c_abs2(g_psi[j]);
    g_mp[j] = 0.0;
  }

  for (int l = 0; l < WG_N; l++) {
    for (int j = 0; j < WG_N; j++) {
      double v = g_w[l * WG_N + j];

      g_mp[l] += v * WG_DX;
      g_wmin = fmin(g_wmin, v);
      g_wmax = fmax(g_wmax, v);
    }
  }

  for (int n = 0; n <= WG_NMAX; n++) {
    g_pn[n] = c_abs2(g_c[n]);
  }

  g_nbar = ps_mean_photon_number(g_c, WG_NMAX);
  g_neg = ps_wigner_negativity(g_w, WG_N, WG_N, WG_DX, WG_DX);

  return 0;
}
