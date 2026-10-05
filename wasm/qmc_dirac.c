/*
 * qmc_dirac.c - WebAssembly front end for the Dirac Klein-paradox and
 * Zitterbewegung demo.
 *
 * All physics comes from the library's physics/dirac_evolve.c: packet
 * construction (dirac_packet_1d), split-operator propagation
 * (dirac_evolve_1d), observables, and the analytic step transmission
 * (dirac_step_transmission). This file only owns one global simulation
 * (N = 2048 points, dx = 0.1, units hbar = m = c = 1) and the scratch buffer
 * JavaScript reads the densities from.
 *
 * Scenes:
 *   qmc_dc_klein(k0, V0)  positive-energy packet at x = -50 hitting a sharp
 *                         step of height V0 at x = 0
 *   qmc_dc_zitter()       wide packet of spinor (1, i)/sqrt(2) at rest, no
 *                         potential
 */

#include "../physics/dirac_evolve.h"
#include "complex.h"
#include "vector.h"
#include <math.h>
#include <stdlib.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { DC_N = 2048 };

#define DC_DX 0.1
#define DC_DT 0.04

static cvector_t *g_u;
static cvector_t *g_l;
static double g_v[DC_N];
static double g_buf[3 * DC_N]; /* |upper|^2, |lower|^2, V */
static double g_t;

static int ensure_alloc(void) {
  if (!g_u) {
    g_u = cvector_alloc(DC_N);
  }
  if (!g_l) {
    g_l = cvector_alloc(DC_N);
  }

  return g_u && g_l;
}

QMC_EXPORT(qmc_dc_n)
int qmc_dc_n(void) { return DC_N; }

QMC_EXPORT(qmc_dc_dx)
double qmc_dc_dx(void) { return DC_DX; }

QMC_EXPORT(qmc_dc_dt)
double qmc_dc_dt(void) { return DC_DT; }

QMC_EXPORT(qmc_dc_time)
double qmc_dc_time(void) { return g_t; }

/* Scratch buffer: [0,N) |upper|^2, [N,2N) |lower|^2, [2N,3N) V. */
QMC_EXPORT(qmc_dc_buffer)
double *qmc_dc_buffer(void) { return g_buf; }

/* Refresh the scratch buffer from the current state. */
QMC_EXPORT(qmc_dc_density)
void qmc_dc_density(void) {
  if (!g_u || !g_l) {
    return;
  }
  for (int i = 0; i < DC_N; i++) {
    g_buf[i] = c_abs2(g_u->data[i]);
    g_buf[DC_N + i] = c_abs2(g_l->data[i]);
    g_buf[2 * DC_N + i] = g_v[i];
  }
}

/* Klein scene. k0 in (0, 3], V0 in [0, 20]. Returns 0 on success. */
QMC_EXPORT(qmc_dc_klein)
int qmc_dc_klein(double k0, double v0) {
  if (!ensure_alloc() || !(k0 > 0.0 && k0 <= 3.0) || !(v0 >= 0.0 && v0 <= 20.0)) {
    return -1;
  }
  for (int i = 0; i < DC_N; i++) {
    g_v[i] = (i > DC_N / 2) ? v0 : 0.0;
  }
  g_t = 0.0;

  return dirac_packet_1d(g_u, g_l, DC_DX, -50.0, k0, 6.0, 1, 1.0, 1.0, 1.0);
}

/* Zitterbewegung scene. Returns 0 on success. */
QMC_EXPORT(qmc_dc_zitter)
int qmc_dc_zitter(void) {
  if (!ensure_alloc()) {
    return -1;
  }
  for (int i = 0; i < DC_N; i++) {
    g_v[i] = 0.0;
  }
  g_t = 0.0;
  if (dirac_packet_1d(g_u, g_l, DC_DX, 0.0, 0.0, 5.0, 0, 1.0, 1.0, 1.0) != 0) {
    return -1;
  }
  for (int i = 0; i < DC_N; i++) {
    complex_t g = g_u->data[i];

    g_u->data[i] = c_scale(g, 1.0 / sqrt(2.0));
    g_l->data[i] = c_scale(c_new(-g.im, g.re), 1.0 / sqrt(2.0));
  }

  return 0;
}

/* Advance n steps of DC_DT. Returns 0 on success. */
QMC_EXPORT(qmc_dc_step)
int qmc_dc_step(int n) {
  if (!g_u || !g_l || n < 1) {
    return -1;
  }
  int rc = dirac_evolve_1d(g_u, g_l, g_v, DC_DX, DC_DT, n, 1.0, 1.0, 1.0);

  if (rc == 0) {
    g_t += DC_DT * n;
  }

  return rc;
}

QMC_EXPORT(qmc_dc_norm)
double qmc_dc_norm(void) { return dirac_norm_1d(g_u, g_l, DC_DX); }

QMC_EXPORT(qmc_dc_position)
double qmc_dc_position(void) { return dirac_position_1d(g_u, g_l, DC_DX); }

/* Probability to the right of the step (x > 0). */
QMC_EXPORT(qmc_dc_right)
double qmc_dc_right(void) {
  if (!g_u || !g_l) {
    return NAN;
  }
  double sum = 0.0;

  for (int i = DC_N / 2 + 1; i < DC_N; i++) {
    sum += c_abs2(g_u->data[i]) + c_abs2(g_l->data[i]);
  }

  return sum * DC_DX;
}

/* Analytic plane-wave transmission for wave number k0 and step V0. */
QMC_EXPORT(qmc_dc_exact)
double qmc_dc_exact(double k0, double v0) {
  return dirac_step_transmission(sqrt(k0 * k0 + 1.0), v0, 1.0, 1.0, 1.0);
}
