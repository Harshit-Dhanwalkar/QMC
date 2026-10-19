/*
 * qmc_corral.c - quantum-corral WebAssembly demo
 *
 * Physics: physics/circular_billiard.c (Bessel-function eigenmodes of the unit
 * disc). A quantum corral is a 2D electron gas confined to a circular region:
 * its eigenstates are Bessel functions J_m(k_mn r) \exp(i m \theta), with
 * energies proportional to k_mn^2. The wavefunction of a launched packet is
 * written as a superposition of these modes; the time evolution multiplies
 * each coefficient by \exp(-i E_mn t), so the field at any later time is a
 * single pass over the retained modes
 */

#include "../physics/circular_billiard.h"
#include <math.h>
#include <stddef.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { QB_N = 80, QB_MAX_MODES = 192 };

// Upper limit on |k| for retained modes
// Number of modes with |k| <= K grows as K^2 (2D density of states), so this
// cutoff is what sets both cost of one-time table build and size of static
// table
static const double QB_KMAX = 26.0;

static cb_mode_t g_modes[QB_MAX_MODES];
static int g_count;
static double g_table[QB_MAX_MODES * QB_N *
                      QB_N]; /* mode value on the grid, 0 outside the disc */
static unsigned char g_inside[QB_N * QB_N];
static double g_coef_re[QB_MAX_MODES], g_coef_im[QB_MAX_MODES];
static double
    g_field[QB_N * QB_N]; /* last field written: mode value or packet density */
static double g_leak, g_energy;
static int g_ready, g_has_packet;

// Map pixel (ix, iy) to its Cartesian centre in [-1, 1]^2. Half-pixel offset
// matters: grid point is centre of pixel, not its corner, which keeps every
// quadrature sum symmetric under reflection
static void grid_point(int ix, int iy, double *x, double *y) {
  *x = -1.0 + 2.0 * (ix + 0.5) / QB_N;
  *y = -1.0 + 2.0 * (iy + 0.5) / QB_N;
}

// Lazy one-time setup: enumerate modes, build inside/outside mask and fill mode
// table. Idempotent
static void init(void) {
  if (g_ready) {
    return;
  }

  g_count = cb_modes(QB_KMAX, g_modes, QB_MAX_MODES);
  if (g_count > QB_MAX_MODES) {
    g_count = QB_MAX_MODES;
  }

  for (int iy = 0; iy < QB_N; iy++) {
    for (int ix = 0; ix < QB_N; ix++) {
      double x;
      double y;

      grid_point(ix, iy, &x, &y);
      g_inside[iy * QB_N + ix] = hypot(x, y) < 1.0;
    }
  }

  for (int i = 0; i < g_count; i++) {
    for (int iy = 0; iy < QB_N; iy++) {
      for (int ix = 0; ix < QB_N; ix++) {
        double x;
        double y;
        double r;

        grid_point(ix, iy, &x, &y);
        r = hypot(x, y);
        g_table[(size_t)i * QB_N * QB_N + iy * QB_N + ix] =
            r < 1.0 ? cb_mode_value(&g_modes[i], r, atan2(y, x)) : 0.0;
      }
    }
  }

  g_ready = 1;
}

// Grid width: QB_N pixels on each side
QMC_EXPORT(qmc_cb_n)
int qmc_cb_n(void) { return QB_N; }

// Number of retained modes with |k| <= QB_KMAX
// Forces lazy setup
QMC_EXPORT(qmc_cb_mode_count)
int qmc_cb_mode_count(void) {
  init();

  return g_count;
}

// last field drawn (a mode or packet density), row-major, y rows from -1 to 1,
// 0 outside disc
// Length QB_N * QB_N
QMC_EXPORT(qmc_cb_field_buffer)
double *qmc_cb_field_buffer(void) { return g_field; }

// Angular quantum number of mode i (sorted by energy), or -1 on bad index
QMC_EXPORT(qmc_cb_mode_m)
int qmc_cb_mode_m(int i) {
  init();

  return i >= 0 && i < g_count ? g_modes[i].m : -1;
}

// Radial quantum number of mode i, or -1 on bad index
QMC_EXPORT(qmc_cb_mode_n)
int qmc_cb_mode_n(int i) {
  init();

  return i >= 0 && i < g_count ? g_modes[i].n : -1;
}

// 1 if mode i uses \sin(m \theta) (odd m), 0 if \cos(m \theta); -1 on bad index
QMC_EXPORT(qmc_cb_mode_sine)
int qmc_cb_mode_sine(int i) {
  init();

  return i >= 0 && i < g_count ? g_modes[i].sine : -1;
}

// Wavenumber k of mode i, or NaN on bad index
QMC_EXPORT(qmc_cb_mode_k)
double qmc_cb_mode_k(int i) {
  init();

  return i >= 0 && i < g_count ? g_modes[i].k : NAN;
}

// Copy mode i into field buffer. Returns 0, or -1 on bad index
QMC_EXPORT(qmc_cb_show_mode)
int qmc_cb_show_mode(int i) {
  init();

  if (i < 0 || i >= g_count) {
    return -1;
  }

  for (int p = 0; p < QB_N * QB_N; p++) {
    g_field[p] = g_table[(size_t)i * QB_N * QB_N + p];
  }

  return 0;
}

/*
 * Launch a Gaussian packet \exp(-|r - r0|^2 / (4 s^2) + i p.r) centred at (x0,
 * y0), with spread s in [0.05, 0.4] and momentum (px, py) with |p| <= 20
 * Packet must start at least 0.05 inside disc boundary (|r0| < 0.95); closer to
 * wall, truncation error from sharp boundary is large
 *
 * Coefficients are overlaps c_i = <\phi_i | \psi> with every mode, evaluated by
 * a Riemann sum on same grid modes live on, then renormalised so that \sum
 * |c_i|^2 = 1 on retained subspace
 *
 * Returns fraction of packet's weight captured by retained modes (1 = exact),
 * or -1 on bad input
 */
QMC_EXPORT(qmc_cb_launch)
double qmc_cb_launch(double x0, double y0, double s, double px, double py) {
  init();
  if (!(hypot(x0, y0) < 0.95) || !(s >= 0.05 && s <= 0.4) ||
      !(hypot(px, py) <= 20.0)) {
    return -1.0;
  }

  double h = 2.0 / QB_N;
  double total = 0.0;
  double kept = 0.0;

  // Norm of Gaussian envelope restricted to disc. This is denominator that
  // normalises packet; it will differ from free-space value 2 \pi s^2 when
  // packet overlaps boundary and that difference is exactly what leak measures
  for (int iy = 0; iy < QB_N; iy++) {
    for (int ix = 0; ix < QB_N; ix++) {
      double x;
      double y;

      grid_point(ix, iy, &x, &y);
      if (g_inside[iy * QB_N + ix]) {
        total +=
            exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (2.0 * s * s)) *
            h * h;
      }
    }
  }

  // Overlap with each mode
  // Integrand is mode value times Gaussian envelope times plane wave; sum over
  // grid is Riemann approximation to <\phi_i | \psi>
  for (int i = 0; i < g_count; i++) {
    double re = 0.0, im = 0.0;
    const double *tab = &g_table[(size_t)i * QB_N * QB_N];

    for (int iy = 0; iy < QB_N; iy++) {
      for (int ix = 0; ix < QB_N; ix++) {
        double x;
        double y;
        double g;
        double ph;

        grid_point(ix, iy, &x, &y);
        if (!g_inside[iy * QB_N + ix]) {
          continue;
        }

        g = exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (4.0 * s * s)) *
            tab[iy * QB_N + ix] * h * h;

        ph = px * x + py * y;
        re += g * cos(ph);
        im += g * sin(ph);
      }
    }

    g_coef_re[i] = re / sqrt(total);
    g_coef_im[i] = im / sqrt(total);
    kept += g_coef_re[i] * g_coef_re[i] + g_coef_im[i] * g_coef_im[i];
  }

  // Mean retained energy <k^2>, then renormalise coefficients to unit weight on
  // retained subspace. 2 are consistent: kept in energy accumulator matches
  // \sqrt(kept) scaling applied to  coefficients in same loop
  g_energy = 0.0;
  for (int i = 0; i < g_count; i++) {
    g_energy += g_modes[i].k * g_modes[i].k *
                (g_coef_re[i] * g_coef_re[i] + g_coef_im[i] * g_coef_im[i]) /
                kept;
    g_coef_re[i] /= sqrt(kept);
    g_coef_im[i] /= sqrt(kept);
  }

  g_leak = kept;
  g_has_packet = 1;

  return kept;
}

// Mean energy <k^2> of launched packet
// NaN before a launch
QMC_EXPORT(qmc_cb_energy)
double qmc_cb_energy(void) { return g_has_packet ? g_energy : NAN; }

// Fill field buffer with |\psi(x, y, t)|^2 of launched packet; returns total
// probability on grid (~1)
QMC_EXPORT(qmc_cb_evolve)
double qmc_cb_evolve(double t) {
  if (!g_has_packet || !isfinite(t)) {
    return NAN;
  }

  double re[QB_N * QB_N];
  double im[QB_N * QB_N];
  double total = 0.0;
  double h = 2.0 / QB_N;

  for (int p = 0; p < QB_N * QB_N; p++) {
    re[p] = im[p] = 0.0;
  }

  // Each mode gets a phase \exp(-i k^2 t). Multiplying coefficient by this
  // phase rotates (re, im) in complex plane; real and imaginary parts of
  // coefficient mix into each other
  for (int i = 0; i < g_count; i++) {
    double a = -g_modes[i].k * g_modes[i].k * t;
    double c = cos(a);
    double s = sin(a);
    double ar = g_coef_re[i] * c - g_coef_im[i] * s;
    double ai = g_coef_re[i] * s + g_coef_im[i] * c;
    const double *tab = &g_table[(size_t)i * QB_N * QB_N];

    for (int p = 0; p < QB_N * QB_N; p++) {
      re[p] += ar * tab[p];
      im[p] += ai * tab[p];
    }
  }

  // Density and total probability on grid. Total is ~1 up to  Riemann error on
  // mode table itself, which vanishes as grid resolution grows
  for (int p = 0; p < QB_N * QB_N; p++) {
    g_field[p] = re[p] * re[p] + im[p] * im[p];
    total += g_field[p] * h * h;
  }

  return total;
}
