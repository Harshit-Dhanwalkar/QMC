/*
 * qmc_wasm.c - thin WebAssembly front end for QMC's 2D split-operator solver
 *
 * Units: hbar = m = 1
 * Grid layout: psi[ix * N + iy], so ix is x coordinate and iy the y
 * coordinate. The same layout is used for potential and density buffers
 */

#include "../core/complex.h"
#include "../core/vector.h"
#include "../physics/soft.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { QMC_MAX_N = 256 };

static int g_n = 0;
static double g_dx = 0.0;
static double g_dt = 0.0;
static double g_box = 0.0;
static cvector_t *g_psi = NULL;
static double *g_v = NULL;
static float *g_density = NULL;
static double *g_mask = NULL;

static double coord(int i) { return (i - g_n / 2) * g_dx; }

static void free_state(void) {
  cvector_free(g_psi);
  free(g_v);
  free(g_density);
  free(g_mask);

  g_psi = NULL;
  g_v = NULL;
  g_density = NULL;
  g_mask = NULL;
  g_n = 0;
}

/* Edge damping so packet leaves box instead of wrapping around (FFT grid is
 * periodic). Applied once per step after soft_evolve_2d() */
static void build_mask(int edge_cells, double strength) {
  for (int ix = 0; ix < g_n; ix++) {
    for (int iy = 0; iy < g_n; iy++) {
      int dx_edge = ix < g_n - 1 - ix ? ix : g_n - 1 - ix;
      int dy_edge = iy < g_n - 1 - iy ? iy : g_n - 1 - iy;
      int d = dx_edge < dy_edge ? dx_edge : dy_edge;
      double f = 1.0;
      if (d < edge_cells) {
        double u = (double)(edge_cells - d) / edge_cells;
        f = 1.0 - strength * u * u;
      }

      g_mask[ix * g_n + iy] = f;
    }
  }
}

/* n must be a power of two (library FFT is radix-2), 8 <= n <= 256.
 * box is side length of simulation square, dt time step
 *
 * Returns 0 on success, -1 on bad arguments or out of memory. */
QMC_EXPORT(qmc_init)
int qmc_init(int n, double box, double dt) {
  if (n < 8 || n > QMC_MAX_N || (n & (n - 1)) != 0 || box <= 0.0 || dt <= 0.0) {
    return -1;
  }

  free_state();

  g_n = n;
  g_box = box;
  g_dx = box / n;
  g_dt = dt;

  size_t cells = (size_t)n * (size_t)n;

  g_psi = cvector_alloc(n * n);
  g_v = calloc(cells, sizeof *g_v);
  g_density = calloc(cells, sizeof *g_density);
  g_mask = calloc(cells, sizeof *g_mask);

  if (!g_psi || !g_v || !g_density || !g_mask) {
    free_state();

    return -1;
  }

  build_mask(n / 12 > 2 ? n / 12 : 2, 0.12);

  return 0;
}

QMC_EXPORT(qmc_n)
int qmc_n(void) { return g_n; }

QMC_EXPORT(qmc_box)
double qmc_box(void) { return g_box; }

/* Pointer to potential V[ix * N + iy] (double). JavaScript may write to it
 * directly (e.g. to paint walls); it is read at every qmc_step() */
QMC_EXPORT(qmc_potential)
double *qmc_potential(void) { return g_v; }

QMC_EXPORT(qmc_clear_potential)
void qmc_clear_potential(void) {
  if (g_v) {
    memset(g_v, 0, (size_t)g_n * (size_t)g_n * sizeof *g_v);
  }
}

/* Fill a rectangle (physical coordinates) with potential value v0 */
QMC_EXPORT(qmc_fill_rect)
void qmc_fill_rect(double x0, double x1, double y0, double y1, double v0) {
  if (!g_v) {
    return;
  }

  for (int ix = 0; ix < g_n; ix++) {
    double x = coord(ix);
    if (x < x0 || x > x1) {
      continue;
    }

    for (int iy = 0; iy < g_n; iy++) {
      double y = coord(iy);
      if (y >= y0 && y <= y1) {
        g_v[ix * g_n + iy] = v0;
      }
    }
  }
}

/* Vertical wall centred at wall_x with given thickness, with `slits` (1 or 2)
 * openings of width slit_w; for 2 slits centre-to-centre separation is sep
 * Wall height (potential) is v0 */
QMC_EXPORT(qmc_build_slits)
void qmc_build_slits(double wall_x, double thickness, int slits, double sep,
                     double slit_w, double v0) {
  if (!g_v) {
    return;
  }

  qmc_clear_potential();
  double half = 0.5 * thickness;
  for (int ix = 0; ix < g_n; ix++) {
    double x = coord(ix);
    if (fabs(x - wall_x) > half) {
      continue;
    }

    for (int iy = 0; iy < g_n; iy++) {
      double y = coord(iy);
      int open = 0;
      if (slits == 1) {
        open = fabs(y) <= 0.5 * slit_w;
      } else {
        open = fabs(fabs(y) - 0.5 * sep) <= 0.5 * slit_w;
      }

      g_v[ix * g_n + iy] = open ? 0.0 : v0;
    }
  }
}

/* Gaussian wavepacket with momentum (kx, ky) and position spread \sigma */
QMC_EXPORT(qmc_wavepacket)
int qmc_wavepacket(double x0, double y0, double kx, double ky, double sigma) {
  if (!g_psi || sigma <= 0.0) {
    return -1;
  }

  for (int ix = 0; ix < g_n; ix++) {
    double x = coord(ix) - x0;
    for (int iy = 0; iy < g_n; iy++) {
      double y = coord(iy) - y0;
      double env = exp(-(x * x + y * y) / (4.0 * sigma * sigma));
      double phase = kx * x + ky * y;
      g_psi->data[ix * g_n + iy] = c_new(env * cos(phase), env * sin(phase));
    }
  }

  double norm = grid_norm(g_psi, g_dx * g_dx);
  if (norm <= 0.0) {
    return -1;
  }

  double s = 1.0 / sqrt(norm);
  for (int i = 0; i < g_n * g_n; i++) {
    g_psi->data[i] = c_scale(g_psi->data[i], s);
  }

  return 0;
}

/* Advance by `steps` time steps
 *
 * Returns 0 on success, -1 on failure
 */
QMC_EXPORT(qmc_step)
int qmc_step(int steps) {
  if (!g_psi) {
    return -1;
  }

  for (int s = 0; s < steps; s++) {
    if (soft_evolve_2d(g_psi, g_v, g_n, g_n, g_dx, g_dx, g_dt, 1, 1.0, 1.0) !=
        0) {
      return -1;
    }

    for (int i = 0; i < g_n * g_n; i++) {
      g_psi->data[i] = c_scale(g_psi->data[i], g_mask[i]);
    }
  }

  return 0;
}

/* Total probability inside box (starts at 1, drops as the packet is absorbed at
 * the edges) */
QMC_EXPORT(qmc_norm)
double qmc_norm(void) { return g_psi ? grid_norm(g_psi, g_dx * g_dx) : 0.0; }

/* Refresh and return probability density |\psi|^2 (float, N*N) */
QMC_EXPORT(qmc_density)
float *qmc_density(void) {
  if (!g_psi) {
    return NULL;
  }

  for (int i = 0; i < g_n * g_n; i++) {
    g_density[i] = (float)c_abs2(g_psi->data[i]);
  }

  return g_density;
}
