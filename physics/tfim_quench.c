#include "tfim_quench.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int valid_input(double J, double h_i, double h_f, double t,
                       int n_sites) {
  return isfinite(J) && isfinite(h_i) && isfinite(h_f) && isfinite(t) &&
         n_sites >= 2 && n_sites <= TFIM_QUENCH_MAX_SITES && n_sites % 2 == 0;
}

// Unit vector of field (0, J \sin(k), h - J \cos(k))
static void field_dir(double J, double h, double k, double n[3]) {
  double y = J * sin(k), z = h - J * cos(k), len = sqrt(y * y + z * z);

  if (len < 1e-300) {
    len = 1.0;
  }

  n[0] = 0.0;
  n[1] = y / len;
  n[2] = z / len;
}

static double quasi_energy(double J, double h, double k) {
  return 2.0 * sqrt(J * J + h * h - 2.0 * J * h * cos(k));
}

static void bloch_at(double J, double h_i, double h_f, double t, double k,
                     double s[3]) {
  double ni[3], nf[3], s0[3];

  field_dir(J, h_i, k, ni);
  field_dir(J, h_f, k, nf);
  for (int c = 0; c < 3; c++) {
    s0[c] = -ni[c];
  }

  double th = 2.0 * quasi_energy(J, h_f, k) * t, co = cos(th), si = sin(th);
  double dot = nf[0] * s0[0] + nf[1] * s0[1] + nf[2] * s0[2];
  double cr[3] = {nf[1] * s0[2] - nf[2] * s0[1], nf[2] * s0[0] - nf[0] * s0[2],
                  nf[0] * s0[1] - nf[1] * s0[0]};

  for (int c = 0; c < 3; c++) {
    s[c] = s0[c] * co + cr[c] * si + nf[c] * dot * (1.0 - co);
  }
}

int tfim_quench_modes(double J, double h_i, double h_f, double t, int n_sites,
                      double *s) {
  if (!s || !valid_input(J, h_i, h_f, t, n_sites)) {
    return -1;
  }

  for (int m = 0; m < n_sites / 2; m++) {
    bloch_at(J, h_i, h_f, t, (2 * m + 1) * M_PI / n_sites, s + 3 * m);
  }

  return 0;
}

double tfim_quench_mx(double J, double h_i, double h_f, double t, int n_sites) {
  if (!valid_input(J, h_i, h_f, t, n_sites)) {
    return NAN;
  }

  double sum = 0.0;
  for (int m = 0; m < n_sites / 2; m++) {
    double s[3];

    bloch_at(J, h_i, h_f, t, (2 * m + 1) * M_PI / n_sites, s);
    sum += s[2];
  }

  return -2.0 * sum / n_sites;
}

double tfim_quench_loschmidt_rate(double J, double h_i, double h_f, double t,
                                  int n_sites) {
  if (!valid_input(J, h_i, h_f, t, n_sites)) {
    return NAN;
  }

  double sum = 0.0;
  for (int m = 0; m < n_sites / 2; m++) {
    double k = (2 * m + 1) * M_PI / n_sites, ni[3], nf[3];

    field_dir(J, h_i, k, ni);
    field_dir(J, h_f, k, nf);
    double c = ni[1] * nf[1] + ni[2] * nf[2],
           sn = sin(quasi_energy(J, h_f, k) * t);
    double lk = 1.0 - (1.0 - c * c) * sn * sn;

    sum += log(lk > 1e-300 ? lk : 1e-300);
  }

  return sum == 0.0 ? 0.0 : -sum / n_sites;
}

/* Pfaffian of a real antisymmetric n x n matrix (n even) by Parlett-Reid
 * elimination with pivoting; destroys a */
static double pfaffian(double *a, int n) {
  double res = 1.0;

  for (int k = 0; k < n - 1; k += 2) {
    int piv = k + 1;

    for (int i = k + 2; i < n; i++) {
      if (fabs(a[i * n + k]) > fabs(a[piv * n + k])) {
        piv = i;
      }
    }

    if (piv != k + 1) {
      for (int j = 0; j < n; j++) {
        double tmp = a[(k + 1) * n + j];

        a[(k + 1) * n + j] = a[piv * n + j];
        a[piv * n + j] = tmp;
      }

      for (int i = 0; i < n; i++) {
        double tmp = a[i * n + k + 1];

        a[i * n + k + 1] = a[i * n + piv];
        a[i * n + piv] = tmp;
      }

      res = -res;
    }

    double p = a[k * n + k + 1];
    if (p == 0.0) {
      return 0.0;
    }

    res *= p;
    for (int i = k + 2; i < n; i++) {
      double tau = a[k * n + i] / p;

      for (int j = k + 2; j < n; j++) {
        a[i * n + j] +=
            tau * a[j * n + k + 1] - a[i * n + k + 1] * (a[k * n + j] / p);
      }
    }
  }

  return res;
}

int tfim_quench_zz(double J, double h_i, double h_f, double t, int n_sites,
                   int r_max, double *zz) {
  if (!zz || !valid_input(J, h_i, h_f, t, n_sites) || r_max < 0 ||
      r_max > TFIM_QUENCH_MAX_R || r_max > n_sites / 2 - 1) {
    return -1;
  }

  int span = 2 * r_max + 3, nm = n_sites / 2;
  double *g = calloc((size_t)span, sizeof *g);
  double *x = calloc((size_t)span, sizeof *x);
  double *rm = calloc((size_t)(4 * r_max * r_max + 1), sizeof *rm);

  if (!g || !x || !rm) {
    free(g);
    free(x);
    free(rm);

    return -2;
  }

  // g[l + off] = G(l), x[l + off] = X(l) for l = -off .. off
  int off = r_max + 1;

  for (int m = 0; m < nm; m++) {
    double k = (2 * m + 1) * M_PI / n_sites, s[3];

    bloch_at(J, h_i, h_f, t, k, s);
    for (int l = -off; l <= off; l++) {
      g[l + off] += (2.0 / n_sites) * (cos(k * l) * -s[2] + sin(k * l) * s[1]);
      x[l + off] += (2.0 / n_sites) * sin(k * l) * s[0];
    }
  }

  zz[0] = 1.0;
  for (int r = 1; r <= r_max; r++) {
    int n = 2 * r;

    // operators B_0 A_1 B_1 A_2 ... B_(r-1) A_r: index 2i = B_i, 2i+1 = A_(i+1)
    for (int a = 0; a < n; a++) {
      rm[a * n + a] = 0.0;
      for (int b = a + 1; b < n; b++) {
        int ia = a / 2 + (a % 2);
        int ib = b / 2 + (b % 2);
        int a_is_a = a % 2;
        int b_is_a = b % 2;
        double v;

        if (a_is_a && b_is_a) {
          v = x[ib - ia + off];
        } else if (!a_is_a && !b_is_a) {
          v = -x[ib - ia + off];
        } else if (!a_is_a) {
          v = g[ib - ia + off];
        } else {
          v = -g[ia - ib + off];
        }

        rm[a * n + b] = v;
        rm[b * n + a] = -v;
      }
    }

    double pf = pfaffian(rm, n);

    zz[r] = (r % 2 ? -pf : pf);
  }

  free(g);
  free(x);
  free(rm);

  return 0;
}

double tfim_quench_critical_time(double J, double h_i, double h_f, int n) {
  if (n < 0 || !isfinite(J) || !isfinite(h_i) || !isfinite(h_f) || J == 0.0 ||
      (h_i - J) * (h_f - J) >= 0.0 || (h_i + h_f) == 0.0) {
    return NAN;
  }

  double c = (J * J + h_i * h_f) / (J * (h_i + h_f));
  if (!(fabs(c) <= 1.0)) {
    return NAN;
  }

  return M_PI * (n + 0.5) / quasi_energy(J, h_f, acos(c));
}
