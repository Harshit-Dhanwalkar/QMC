/*
 * qmc_doublewell.c - WebAssembly front end for double-well tunnelling demo
 *
 * Physics: physics/double_well.c (eigenstates of V = B (1 - x^2 / a^2)^2 on a
 * grid, \hbar = m = 1). Grid has DW_N points on [-DW_L, DW_L]; lowest DW_STATES
 * eigenstates are kept for packet dynamics
 */

#include "../physics/double_well.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { DW_N = 500, DW_STATES = 8, DW_CURVE = 61 };

static const double DW_L = 5.0;
static double g_v[DW_N];
static double g_x[DW_N];
static double g_psi[DW_STATES * DW_N];
static double g_e[DW_STATES];
static double g_c[DW_STATES];
static double g_dens[DW_N];
static double g_re[DW_N];
static double g_im[DW_N];
static double g_split[DW_CURVE];
static double g_wkb[DW_CURVE];
static double g_b = 8.0;
static double g_a = 1.5;
static double g_mean_x;
static int g_ready;
static int g_has_packet;

QMC_EXPORT(qmc_dw_n)
int qmc_dw_n(void) { return DW_N; }
QMC_EXPORT(qmc_dw_states)
int qmc_dw_states(void) { return DW_STATES; }
QMC_EXPORT(qmc_dw_curve_len)
int qmc_dw_curve_len(void) { return DW_CURVE; }
QMC_EXPORT(qmc_dw_half_width)
double qmc_dw_half_width(void) { return DW_L; }

// Grid positions, potential, density of evolved packet (qmc_dw_evolve),
// eigenstates (state k at k * n) and energies
QMC_EXPORT(qmc_dw_x_buffer)
double *qmc_dw_x_buffer(void) { return g_x; }
QMC_EXPORT(qmc_dw_v_buffer)
double *qmc_dw_v_buffer(void) { return g_v; }
QMC_EXPORT(qmc_dw_density_buffer)
double *qmc_dw_density_buffer(void) { return g_dens; }
QMC_EXPORT(qmc_dw_psi_buffer)
double *qmc_dw_psi_buffer(void) { return g_psi; }
QMC_EXPORT(qmc_dw_e_buffer)
double *qmc_dw_e_buffer(void) { return g_e; }

// Splitting \Delta and its WKB estimate versus barrier B = bmax (i + 1) /
// curve_len at current well separation
QMC_EXPORT(qmc_dw_split_buffer)
double *qmc_dw_split_buffer(void) { return g_split; }
QMC_EXPORT(qmc_dw_wkb_buffer)
double *qmc_dw_wkb_buffer(void) { return g_wkb; }

// Barrier height in [1, 60], well centre a in [0.5, 3]. Solves eigenproblem.
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_dw_set)
int qmc_dw_set(double barrier, double a) {
  if (!(barrier >= 1.0 && barrier <= 60.0 && a >= 0.5 && a <= 3.0)) {
    return -1;
  }

  double dx = dw_spacing(DW_N, DW_L);

  if (dw_potential(barrier, a, DW_L, DW_N, g_v) != 0 ||
      dw_eigen(g_v, DW_N, dx, DW_STATES, g_e, g_psi) != 0) {
    return -1;
  }

  for (int i = 0; i < DW_N; i++) {
    g_x[i] = -DW_L + (i + 1) * dx;
  }

  g_b = barrier;
  g_a = a;
  g_ready = 1;
  g_has_packet = 0;

  return 0;
}

QMC_EXPORT(qmc_dw_splitting)
double qmc_dw_splitting(void) { return g_ready ? g_e[1] - g_e[0] : NAN; }

// WKB estimate of splitting for barrier b and separation a: (\omega / \pi)
// \exp(-\int \sqrt(2 (V - E0)) dx)
static double wkb(double b, double a) {
  double w = sqrt(8.0 * b / (a * a));
  double e0 = 0.5 * w;

  if (e0 >= b) {
    return NAN;
  }

  double xt = a * sqrt(1.0 - sqrt(e0 / b));
  double s = 0.0;
  const int m = 600;

  for (int i = 0; i < m; i++) {
    double x = -xt + 2.0 * xt * (i + 0.5) / m;
    double u = 1.0 - x * x / (a * a);

    s += sqrt(fmax(0.0, 2.0 * (b * u * u - e0))) * 2.0 * xt / m;
  }

  return w / 3.141592653589793 * exp(-s);
}

QMC_EXPORT(qmc_dw_wkb)
double qmc_dw_wkb(void) { return wkb(g_b, g_a); }

// Splitting and WKB curves over barrier heights (0, bmax] at current separation
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_dw_split_curve)
int qmc_dw_split_curve(double bmax) {
  static double v[DW_N], psi[2 * DW_N];
  double e[2], dx = dw_spacing(DW_N, DW_L);

  if (!(bmax >= 2.0 && bmax <= 60.0)) {
    return -1;
  }

  for (int i = 0; i < DW_CURVE; i++) {
    double b = bmax * (i + 1) / DW_CURVE;

    if (dw_potential(b, g_a, DW_L, DW_N, v) != 0 ||
        dw_eigen(v, DW_N, dx, 2, e, psi) != 0) {
      return -1;
    }

    g_split[i] = e[1] - e[0];
    g_wkb[i] = wkb(b, g_a);
  }

  return 0;
}

// Launch a Gaussian of width sigma (hbar = m = 1) at `center`
// Returns weight captured by kept states, or -1
QMC_EXPORT(qmc_dw_launch)
double qmc_dw_launch(double center, double sigma) {
  double cap;

  if (!g_ready || !(sigma >= 0.1 && sigma <= 1.5) || !(fabs(center) <= 3.5) ||
      dw_project(g_psi, DW_N, dw_spacing(DW_N, DW_L), DW_STATES, DW_L, center,
                 sigma, g_c, &cap) != 0) {
    return -1.0;
  }

  g_has_packet = 1;

  return cap;
}

// Evolve launched packet to time t; fills density buffer (normalised to its
// captured weight) and returns P(x < 0)
QMC_EXPORT(qmc_dw_evolve)
double qmc_dw_evolve(double t) {
  double dx = dw_spacing(DW_N, DW_L);

  if (!g_has_packet ||
      dw_evolve(g_psi, g_e, g_c, DW_N, DW_STATES, t, g_re, g_im) != 0) {
    return NAN;
  }

  for (int i = 0; i < DW_N; i++) {
    g_dens[i] = g_re[i] * g_re[i] + g_im[i] * g_im[i];
  }

  g_mean_x = dw_mean_x(g_re, g_im, DW_N, dx, DW_L);

  return dw_prob_left(g_re, g_im, DW_N, dx);
}

QMC_EXPORT(qmc_dw_mean_x)
double qmc_dw_mean_x(void) { return g_has_packet ? g_mean_x : NAN; }
