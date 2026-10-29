/*
 * Zeeman effect of a p electron at any field: weak-field (anomalous Zeeman) to
 * Paschen-Back, with line intensities
 *
 * H = \xi L.S + b (L_z + 2 S_z) for l = 1, s = 1/2
 */

#include "zeeman_full.h"
#include <math.h>

// Basis ordering: index (m_l + 1) * 2 + (1 if spin down else 0)
static int idx(int ml, int down) { return (ml + 1) * 2 + down; }

static void clear(zf_state_t *s) {
  for (int i = 0; i < 6; i++) {
    s->amp[i] = 0.0;
  }
}

// Eigen-decomposition of [[a, c], [c, d]]: upper eigenvalue/vector (\cos, \sin)
// and lower (-\sin, \cos)
static void block(double a, double c, double d, double *lo, double *hi,
                  double *cs, double *sn) {
  double mean = 0.5 * (a + d);
  double diff = 0.5 * (a - d);
  double r = sqrt(diff * diff + c * c);
  double theta = 0.5 * atan2(2.0 * c, a - d);

  *hi = mean + r;
  *lo = mean - r;
  *cs = cos(theta);
  *sn = sin(theta);
}

int zf_solve(double xi, double b, zf_state_t s[6]) {
  if (!s || !isfinite(xi) || !isfinite(b)) {
    return -1;
  }

  for (int k = 0; k < 6; k++) {
    clear(&s[k]);
  }

  // m_j = +3/2: |1, up>
  s[0].energy = 0.5 * xi + 2.0 * b;
  s[0].two_mj = 3;
  s[0].j_branch = 3;
  s[0].amp[idx(1, 0)] = 1.0;

  // m_j = -3/2: |-1, dn>
  s[3].energy = 0.5 * xi - 2.0 * b;
  s[3].two_mj = -3;
  s[3].j_branch = 3;
  s[3].amp[idx(-1, 1)] = 1.0;

  // m_j = +1/2 block in (|1,dn>, |0,up>)
  {
    double lo;
    double hi;
    double cs;
    double sn;
    double c = xi / sqrt(2.0);

    // Matrix is [[a, c],[c, d]] with a = -xi/2 (|1,dn>), d = b (|0,up>);
    // upper vector (cs, sn) over (|1,dn>, |0,up>)
    block(-0.5 * xi, c, b, &lo, &hi, &cs, &sn);

    s[1].energy = hi; // j = 3/2 branch: upper level at b = 0
    s[1].two_mj = 1;
    s[1].j_branch = 3;
    s[1].amp[idx(1, 1)] = cs;
    s[1].amp[idx(0, 0)] = sn;

    s[4].energy = lo;
    s[4].two_mj = 1;
    s[4].j_branch = 1;
    s[4].amp[idx(1, 1)] = -sn;
    s[4].amp[idx(0, 0)] = cs;
  }

  // m_j = -1/2 block in (|0,dn>, |-1,up>)
  {
    double lo;
    double hi;
    double cs;
    double sn;
    double c = xi / sqrt(2.0);

    block(-b, c, -0.5 * xi, &lo, &hi, &cs, &sn);

    s[2].energy = hi;
    s[2].two_mj = -1;
    s[2].j_branch = 3;
    s[2].amp[idx(0, 1)] = cs;
    s[2].amp[idx(-1, 0)] = sn;

    s[5].energy = lo;
    s[5].two_mj = -1;
    s[5].j_branch = 1;
    s[5].amp[idx(0, 1)] = -sn;
    s[5].amp[idx(-1, 0)] = cs;
  }

  // At very large |b| "upper"/"lower" labels swap relative to b = 0 j labels
  // (levels do not cross within a block, so branch that is j = 3/2 at b = 0
  // remains upper one). Nothing to do: labels follow branches.

  return 0;
}

double zf_s_energy(double b, int ms_down) { return ms_down ? -b : b; }

double zf_line_intensity(const zf_state_t *s, int ms_down, int q) {
  if (!s || q < -1 || q > 1 || (ms_down != 0 && ms_down != 1)) {
    return NAN;
  }

  double a = s->amp[idx(q, ms_down)];

  return a * a;
}

double zf_magnetic_moment(const zf_state_t *s) {
  double m = 0.0;

  if (!s) {
    return NAN;
  }

  for (int ml = -1; ml <= 1; ml++) {
    for (int down = 0; down < 2; down++) {
      double a = s->amp[idx(ml, down)];

      m += a * a * (ml + 2.0 * (down ? -0.5 : 0.5));
    }
  }

  return m;
}
