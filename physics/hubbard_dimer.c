/*
 * Exact diagonalisation of two-site Hubbard model
 *
 * Two-site lattice ("Hubbard dimer") with hopping t between sites and an
 * on-site Coulomb repulsion U. Hamiltonian is:
 *
 *   H = -t \sum_s (c^\dagger_{1s} c_{2s} + h.c.)
 *       + U \sum_i n_{i\uparrow} n_{i\downarrow}
 *
 * 16-dimensional Fock space is spanned by 4-bit masks encoding which of 4
 * spin-orbitals {1up, 1dn, 2up, 2dn} are occupied. Hopping terms pick up
 * fermionic Jordan-Wigner signs computed via popcount, so basis is fully
 * antisymmetric
 */

#include "hubbard_dimer.h"
#include <math.h>
#include <string.h>

static int popcount(int x) {
  int c = 0;

  while (x) {
    c += x & 1;
    x >>= 1;
  }

  return c;
}

// out += c+_i c_j |in> (real amplitude vectors over masks)
static void hop(const double *in, double *out, int i, int j) {
  for (int m = 0; m < HD_DIM; m++) {
    if (in[m] == 0.0 || !((m >> j) & 1)) {
      continue;
    }

    int sign = (popcount(m & ((1 << j) - 1)) & 1) ? -1 : 1;
    int m1 = m ^ (1 << j);

    if ((m1 >> i) & 1) {
      continue;
    }

    if (popcount(m1 & ((1 << i) - 1)) & 1) {
      sign = -sign;
    }

    out[m1 | (1 << i)] += sign * in[m];
  }
}

// out = H in
static void apply_h(double t, double u, const double *in, double *out) {
  memset(out, 0, sizeof(double) * HD_DIM);
  for (int s = 0; s < 2; s++) {
    hop(in, out, s, 2 + s); // c+_{1s} c_{2s}
    hop(in, out, 2 + s, s);
  }

  for (int p = 0; p < HD_DIM; p++) {
    out[p] *= -t;
  }

  for (int m = 0; m < HD_DIM; m++) {
    out[m] += u * (((m & 3) == 3) + (((m >> 2) & 3) == 3)) * in[m];
  }
}

// Cyclic Jacobi for a symmetric matrix: a becomes diagonal, columns of vec
// eigenvectors
static void jacobi(double a[HD_DIM][HD_DIM], double vec[HD_DIM][HD_DIM]) {
  for (int i = 0; i < HD_DIM; i++) {
    for (int j = 0; j < HD_DIM; j++) {
      vec[i][j] = i == j;
    }
  }

  for (int sweep = 0; sweep < 60; sweep++) {
    double off = 0.0;

    for (int p = 0; p < HD_DIM; p++) {
      for (int q = p + 1; q < HD_DIM; q++) {
        off += a[p][q] * a[p][q];
      }
    }

    if (off < 1e-30) {
      break;
    }

    for (int p = 0; p < HD_DIM; p++) {
      for (int q = p + 1; q < HD_DIM; q++) {
        if (fabs(a[p][q]) < 1e-300) {
          continue;
        }

        double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        double tt = (theta >= 0 ? 1.0 : -1.0) /
                    (fabs(theta) + sqrt(theta * theta + 1.0));
        double c = 1.0 / sqrt(tt * tt + 1.0);
        double s = tt * c;

        for (int k = 0; k < HD_DIM; k++) {
          double akp = a[k][p];
          double akq = a[k][q];

          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }

        for (int k = 0; k < HD_DIM; k++) {
          double apk = a[p][k];
          double aqk = a[q][k];

          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }

        for (int k = 0; k < HD_DIM; k++) {
          double vkp = vec[k][p];
          double vkq = vec[k][q];

          vec[k][p] = c * vkp - s * vkq;
          vec[k][q] = s * vkp + c * vkq;
        }
      }
    }
  }
}

int hd_solve(double t, double u, hd_solution_t *s) {
  if (!s || !isfinite(t) || !isfinite(u)) {
    return -1;
  }

  double a[HD_DIM][HD_DIM];
  double vec[HD_DIM][HD_DIM];
  double e[HD_DIM];
  double unit[HD_DIM];
  double col[HD_DIM];

  for (int m = 0; m < HD_DIM; m++) {
    memset(unit, 0, sizeof unit);
    unit[m] = 1.0;
    apply_h(t, u, unit, col);

    for (int r = 0; r < HD_DIM; r++) {
      a[r][m] = col[r];
    }
  }

  jacobi(a, vec);
  for (int k = 0; k < HD_DIM; k++) {
    e[k] = a[k][k];
  }

  // sort ascending, ties by particle number then 2 S_z so degenerate multiplets
  // come out in a fixed order
  int order[HD_DIM];
  int key_n[HD_DIM];
  int key_s[HD_DIM];

  for (int k = 0; k < HD_DIM; k++) {
    double nn = 0.0;
    double ss = 0.0;

    for (int m = 0; m < HD_DIM; m++) {
      double w = vec[m][k] * vec[m][k];

      nn += w * popcount(m);
      ss += w * (popcount(m & 5) - popcount(m & 10));
    }

    key_n[k] = (int)lround(nn);
    key_s[k] = (int)lround(ss);
    order[k] = k;
  }

  for (int i = 1; i < HD_DIM; i++) {
    int x = order[i];
    int j = i - 1;

    while (j >= 0) {
      int y = order[j];
      int later = e[y] > e[x] + 1e-9 ||
                  (fabs(e[y] - e[x]) <= 1e-9 &&
                   (key_n[y] > key_n[x] ||
                    (key_n[y] == key_n[x] && key_s[y] > key_s[x])));

      if (!later) {
        break;
      }

      order[j + 1] = y;
      j--;
    }

    order[j + 1] = x;
  }

  s->t = t;
  s->u = u;
  for (int k = 0; k < HD_DIM; k++) {
    int o = order[k];

    s->e[k] = e[o];
    s->n[k] = key_n[o];
    s->twosz[k] = key_s[o];
    for (int m = 0; m < HD_DIM; m++) {
      s->v[k][m] = vec[m][o];
    }
  }

  return 0;
}

void hd_n2_states(const hd_solution_t *s, int idx[6]) {
  int c = 0;

  for (int k = 0; k < HD_DIM && c < 6; k++) {
    if (s->n[k] == 2) {
      idx[c++] = k;
    }
  }
}

void hd_basis_state(int mask, double *re, double *im) {
  for (int m = 0; m < HD_DIM; m++) {
    re[m] = im[m] = 0.0;
  }

  if (mask >= 0 && mask < HD_DIM) {
    re[mask] = 1.0;
  }
}

void hd_eigenstate(const hd_solution_t *s, int k, double *re, double *im) {
  for (int m = 0; m < HD_DIM; m++) {
    re[m] = (k >= 0 && k < HD_DIM) ? s->v[k][m] : 0.0;
    im[m] = 0.0;
  }
}

int hd_evolve(const hd_solution_t *s, const double *re0, const double *im0,
              double time, double *re, double *im) {
  if (!isfinite(time)) {
    return -1;
  }

  double tr[HD_DIM] = {0};
  double ti[HD_DIM] = {0};

  for (int k = 0; k < HD_DIM; k++) {
    // <k|\psi0> (eigenvectors real)
    double cr = 0.0;
    double ci = 0.0;

    for (int m = 0; m < HD_DIM; m++) {
      cr += s->v[k][m] * re0[m];
      ci += s->v[k][m] * im0[m];
    }

    double c = cos(s->e[k] * time);
    double sn = sin(s->e[k] * time);

    // \exp(-iEt (cr + i ci))
    double ar = cr * c + ci * sn;
    double ai = ci * c - cr * sn;

    for (int m = 0; m < HD_DIM; m++) {
      tr[m] += ar * s->v[k][m];
      ti[m] += ai * s->v[k][m];
    }
  }

  for (int m = 0; m < HD_DIM; m++) {
    re[m] = tr[m];
    im[m] = ti[m];
  }

  return 0;
}

double hd_expect_energy(const hd_solution_t *s, const double *re,
                        const double *im) {
  double hr[HD_DIM];
  double hi[HD_DIM];
  double v = 0.0;

  apply_h(s->t, s->u, re, hr);
  apply_h(s->t, s->u, im, hi);
  for (int m = 0; m < HD_DIM; m++) {
    v += re[m] * hr[m] + im[m] * hi[m];
  }

  return v;
}

double hd_expect_double_occ(const double *re, const double *im) {
  double v = 0.0;

  for (int m = 0; m < HD_DIM; m++) {
    v += (re[m] * re[m] + im[m] * im[m]) *
         (((m & 3) == 3) + (((m >> 2) & 3) == 3));
  }

  return v;
}

double hd_expect_sz(const double *re, const double *im, int site) {
  double v = 0.0;
  int up = site ? 4 : 1;
  int dn = site ? 8 : 2;

  for (int m = 0; m < HD_DIM; m++) {
    v += (re[m] * re[m] + im[m] * im[m]) * 0.5 *
         (((m & up) != 0) - ((m & dn) != 0));
  }

  return v;
}

// <S1 . S2> = <Sz1 Sz2> + (1/2) <S1+ S2- + S1- S2+>, S+_i = c+_{i up} c_{i dn}
double hd_expect_ss(const double *re, const double *im) {
  double v = 0.0;
  double ar[HD_DIM] = {0};
  double ai[HD_DIM] = {0};

  for (int m = 0; m < HD_DIM; m++) {
    double sz1 = 0.5 * (((m & 1) != 0) - ((m & 2) != 0));
    double sz2 = 0.5 * (((m & 4) != 0) - ((m & 8) != 0));

    v += (re[m] * re[m] + im[m] * im[m]) * sz1 * sz2;
  }

  // S1+ S2- = c+_{1up} c_{1dn} c+_{2dn} c_{2up}; apply right-hand pair first,
  // then left
  double tr[HD_DIM] = {0};
  double ti[HD_DIM] = {0};
  double ur[HD_DIM] = {0};
  double ui[HD_DIM] = {0};

  hop(re, tr, 3, 2); /* c+_{2dn} c_{2up} */
  hop(im, ti, 3, 2);
  hop(tr, ar, 0, 1); /* c+_{1up} c_{1dn} */
  hop(ti, ai, 0, 1);
  hop(re, ur, 2, 3); /* S1- S2+ = c+_{1dn} c_{1up} c+_{2up} c_{2dn} */
  hop(im, ui, 2, 3);

  double br[HD_DIM] = {0};
  double bi[HD_DIM] = {0};

  hop(ur, br, 1, 0);
  hop(ui, bi, 1, 0);
  for (int m = 0; m < HD_DIM; m++) {
    v += 0.5 * (re[m] * (ar[m] + br[m]) + im[m] * (ai[m] + bi[m]));
  }

  return v;
}

void hd_site_probs(const double *re, const double *im, int site, double p[4]) {
  p[0] = p[1] = p[2] = p[3] = 0.0;
  for (int m = 0; m < HD_DIM; m++) {
    int bits = (m >> (2 * site)) & 3; // bit0 up, bit1 down

    p[bits == 0   ? 0
      : bits == 1 ? 1
      : bits == 2 ? 2
                  : 3] += re[m] * re[m] + im[m] * im[m];
  }
}

double hd_site_entropy(const double *re, const double *im, int site) {
  double p[4];
  double v = 0.0;

  hd_site_probs(re, im, site, p);
  for (int i = 0; i < 4; i++) {
    if (p[i] > 1e-15) {
      v -= p[i] * log(p[i]);
    }
  }

  return v;
}

int hd_thermal_n2(const hd_solution_t *s, double temp, hd_thermal_t *out) {
  if (!s || !out || !(temp > 0.0) || !isfinite(temp)) {
    return -1;
  }

  int idx[6];
  double zero[HD_DIM] = {0};
  double w[6];
  double z = 0.0;
  double e0;
  double re[HD_DIM];

  hd_n2_states(s, idx);
  e0 = s->e[idx[0]];
  out->energy = out->double_occ = out->ss = out->entropy = 0.0;
  for (int i = 0; i < 6; i++) {
    w[i] = exp(-(s->e[idx[i]] - e0) / temp);
    z += w[i];
  }

  for (int i = 0; i < 6; i++) {
    double p = w[i] / z;

    hd_eigenstate(s, idx[i], re, zero);
    out->energy += p * s->e[idx[i]];
    out->double_occ += p * hd_expect_double_occ(re, zero);
    out->ss += p * hd_expect_ss(re, zero);
    out->entropy -= p * log(p);
  }

  return 0;
}
