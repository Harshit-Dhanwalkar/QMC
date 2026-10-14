/*
 * qmc_carpet.c - WebAssembly front end for the quantum-carpet demo.
 *
 * Physics: physics/box_revival.c (Gaussian packet in an infinite square well,
 * exact revival structure). This file owns the packet, a 192-point position
 * grid, the space-time density carpet over one revival time, a density slice
 * at a chosen time and the return-probability curve.
 */

#include "../physics/box_revival.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { QC_NX = 192, QC_NT = 384, QC_MODES = 160, QC_CURVE = 513 };

static complex_t g_c[QC_MODES];
static double g_x[QC_NX];
static double g_carpet[QC_NT * QC_NX];
static double g_slice[QC_NX];
static double g_curve[QC_CURVE];
static double g_leak;
static int g_ready;

QMC_EXPORT(qmc_qc_nx)
int qmc_qc_nx(void) { return QC_NX; }

QMC_EXPORT(qmc_qc_nt)
int qmc_qc_nt(void) { return QC_NT; }

QMC_EXPORT(qmc_qc_curve_len)
int qmc_qc_curve_len(void) { return QC_CURVE; }

/* carpet[it * NX + ix] = |psi(x_ix, t_it)|^2 with t_it = T it / (NT - 1), T = 2 pi. */
QMC_EXPORT(qmc_qc_carpet_buffer)
double *qmc_qc_carpet_buffer(void) { return g_carpet; }

QMC_EXPORT(qmc_qc_slice_buffer)
double *qmc_qc_slice_buffer(void) { return g_slice; }

/* Return probability |<psi_0|psi(t)>|^2 at t = T i / (CURVE - 1). */
QMC_EXPORT(qmc_qc_curve_buffer)
double *qmc_qc_curve_buffer(void) { return g_curve; }

/* Weight of the packet inside the retained modes (1 = exact). */
QMC_EXPORT(qmc_qc_weight)
double qmc_qc_weight(void) { return g_leak; }

/* Mean energy in units of the ground-state energy, and the dominant mode number. */
QMC_EXPORT(qmc_qc_energy)
double qmc_qc_energy(void) { return g_ready ? box_mean_energy(g_c, QC_MODES) : NAN; }

/*
 * Set the packet: centre x0 in (0, 1), spread sigma in [0.02, 0.2], central mode number n0 in [-60, 60].
 * Recomputes the whole carpet and return curve. Returns 0 on success, -1 on bad input.
 */
QMC_EXPORT(qmc_qc_set)
int qmc_qc_set(double x0, double sigma, double n0) {
  if (!(sigma >= 0.02 && sigma <= 0.2) || !(fabs(n0) <= 60.0)) {
    return -1;
  }
  g_leak = box_gaussian_coefficients(g_c, QC_MODES, x0, sigma, n0);
  if (g_leak < 0.0) {
    g_ready = 0;

    return -1;
  }
  g_ready = 1;
  for (int j = 0; j < QC_NX; j++) {
    g_x[j] = (j + 0.5) / QC_NX;
  }
  if (box_density_map(g_c, QC_MODES, g_x, QC_NX, BOX_REVIVAL_TIME, QC_NT, g_carpet) != 0) {
    g_ready = 0;

    return -1;
  }
  for (int i = 0; i < QC_CURVE; i++) {
    g_curve[i] = c_abs2(box_overlap(g_c, QC_MODES, BOX_REVIVAL_TIME * i / (QC_CURVE - 1)));
  }

  return 0;
}

/* Fill the slice buffer with |psi(x, t)|^2 for t = frac * T. Returns the return probability, NaN on bad input. */
QMC_EXPORT(qmc_qc_slice)
double qmc_qc_slice(double frac) {
  complex_t psi[QC_NX];

  if (!g_ready || !isfinite(frac) || box_wavefunction(g_c, QC_MODES, BOX_REVIVAL_TIME * frac, g_x, QC_NX, psi) != 0) {
    return NAN;
  }
  for (int j = 0; j < QC_NX; j++) {
    g_slice[j] = c_abs2(psi[j]);
  }

  return c_abs2(box_overlap(g_c, QC_MODES, BOX_REVIVAL_TIME * frac));
}
