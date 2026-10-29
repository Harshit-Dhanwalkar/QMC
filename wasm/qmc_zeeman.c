/*
 * qmc_zeeman.c - WebAssembly front end for full-field Zeeman demo (2S and 2P of
 * an alkali atom)
 *
 * Physics: physics/zeeman_full.c (H = xi L.S + b (L_z + 2 S_z), closed-form
 * blocks) Energies are in unit of \xi and b
 */

#include "../physics/zeeman_full.h"
#include <math.h>
#include <stddef.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) // native builds (tests) need no export attribute
#endif

enum { ZM_MAX_POINTS = 512, ZM_LINES = 12 };

static double g_curves[8 * ZM_MAX_POINTS]; // level k at k * npts: 6 P states
                                           // then 2S up, 2S down
static double
    g_lines[ZM_LINES * 5]; // energy, intensity, q, P-state index, ms_down
static double g_mj[8];

QMC_EXPORT(qmc_zm_max_points)
int qmc_zm_max_points(void) { return ZM_MAX_POINTS; }
QMC_EXPORT(qmc_zm_curves_buffer)
double *qmc_zm_curves_buffer(void) { return g_curves; }
QMC_EXPORT(qmc_zm_lines_buffer)
double *qmc_zm_lines_buffer(void) { return g_lines; }
// 2 m_j of eight curves (2S entries are \pm 1)
QMC_EXPORT(qmc_zm_mj_buffer)
double *qmc_zm_mj_buffer(void) { return g_mj; }

// Energies of six 2P states and two 2S states for b = bmax i / (npts - 1)
// Returns 0, or -1 on bad input
QMC_EXPORT(qmc_zm_curves)
int qmc_zm_curves(double xi, double bmax, int npts) {
  if (!(xi > 0.0 && xi < 1e6 && bmax > 0.0 && bmax < 1e7) || npts < 2 ||
      npts > ZM_MAX_POINTS) {
    return -1;
  }

  zf_state_t s[6];

  for (int i = 0; i < npts; i++) {
    double b = bmax * i / (npts - 1);

    if (zf_solve(xi, b, s) != 0) {
      return -1;
    }

    for (int k = 0; k < 6; k++) {
      g_curves[k * npts + i] = s[k].energy;
    }
    g_curves[6 * npts + i] = zf_s_energy(b, 0);
    g_curves[7 * npts + i] = zf_s_energy(b, 1);
  }

  zf_solve(xi, 0.0, s);
  for (int k = 0; k < 6; k++) {
    g_mj[k] = s[k].two_mj;
  }
  g_mj[6] = 1;
  g_mj[7] = -1;

  return 0;
}

// 2S -> 2P lines at field b: row j = 6 * ms_down + k is [E_P(k) - E_S(ms),
// intensity, q, k, ms_down] with q = m_l of P state (-1, 0, 1 for \sigma-, \pi,
// \sigma+)
// Returns number of rows with intensity > 1e-12, or -1 on bad input; rows are
// not sorted, and for each (spin, q) intensities over P states sum to
QMC_EXPORT(qmc_zm_lines)
int qmc_zm_lines(double xi, double b) {
  zf_state_t s[6];

  if (!(xi > 0.0 && xi < 1e6 && fabs(b) < 1e7) || zf_solve(xi, b, s) != 0) {
    return -1;
  }

  int count = 0;
  for (int ms = 0; ms < 2; ms++) {
    for (int k = 0; k < 6; k++) {
      int q = (s[k].two_mj - (ms ? -1 : 1)) / 2; /* m_l = m_j - m_s */
      double inten =
          (q >= -1 && q <= 1) ? zf_line_intensity(&s[k], ms, q) : 0.0;
      double *row = &g_lines[((ptrdiff_t)6 * ms + k) * 5];

      row[0] = s[k].energy - zf_s_energy(b, ms);
      row[1] = inten;
      row[2] = q;
      row[3] = k;
      row[4] = ms;
      if (inten > 1e-12) {
        count++;
      }
    }
  }

  return count;
}

// <L_z + 2 S_z> of P state k at field b (slope of its curve), or NaN for bad
// input
QMC_EXPORT(qmc_zm_moment)
double qmc_zm_moment(double xi, double b, int k) {
  zf_state_t s[6];

  if (k < 0 || k > 5 || !(xi > 0.0) || zf_solve(xi, b, s) != 0) {
    return NAN;
  }

  return zf_magnetic_moment(&s[k]);
}
