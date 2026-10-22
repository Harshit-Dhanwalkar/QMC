/*
 * Implementation of Aharonov-Bohm ring solver
 *
 * Gauge convention (flux phase on one bond), analytic persistent-current
 * formula and two-terminal transmission setup
 */

#include "ab_ring.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { AR_MAX = 64 };

double ar_level(int n_sites, int n, double flux, double t) {
  if (n_sites < 1 || !isfinite(flux) || !isfinite(t)) {
    return NAN;
  }

  return -2.0 * t * cos(2.0 * M_PI * (n - flux) / n_sites);
}

int ar_levels_sorted(int n_sites, double flux, double t, double *e) {
  if (n_sites < 1 || n_sites > 4096 || !e || !isfinite(flux) || !isfinite(t)) {
    return -1;
  }

  for (int n = 0; n < n_sites; n++) {
    e[n] = ar_level(n_sites, n, flux, t);
  }

  for (int i = 1; i < n_sites; i++) {
    double x = e[i];
    int j = i - 1;

    while (j >= 0 && e[j] > x) {
      e[j + 1] = e[j];
      j--;
    }

    e[j + 1] = x;
  }

  return 0;
}

// Indices of ne lowest levels, ties (within 1e-12) broken by small distance of
// n
static int lowest(int n_sites, int ne, double flux, double t, int *pick) {
  if (n_sites < 1 || n_sites > 4096 || ne < 0 || ne > n_sites) {
    return -1;
  }

  int used[4096] = {0};

  for (int k = 0; k < ne; k++) {
    int best = -1;
    double eb = 0.0;

    for (int n = 0; n < n_sites; n++) {
      double e = ar_level(n_sites, n, flux, t);

      if (!used[n] && (best < 0 || e < eb - 1e-12)) {
        best = n;
        eb = e;
      }
    }

    // If all sites are used, best remains -1.
    if (best < 0) {
      return -1;
    }

    used[best] = 1;
    pick[k] = best;
  }

  return 0;
}

double ar_energy(int n_sites, int ne, double flux, double t) {
  if (n_sites < 1 || n_sites > 4096 || ne < 0 || ne > n_sites ||
      !isfinite(flux) || !isfinite(t)) {
    return NAN;
  }

  int pick[4096];
  double s = 0.0;

  if (lowest(n_sites, ne, flux, t, pick) != 0) {
    return NAN;
  }

  for (int k = 0; k < ne; k++) {
    s += ar_level(n_sites, pick[k], flux, t);
  }

  return s;
}

double ar_current(int n_sites, int ne, double flux, double t) {
  if (n_sites < 1 || n_sites > 4096 || ne < 0 || ne > n_sites ||
      !isfinite(flux) || !isfinite(t)) {
    return NAN;
  }

  int pick[4096];
  double s = 0.0;

  lowest(n_sites, ne, flux, t, pick);
  for (int k = 0; k < ne; k++) {
    s +=
        4.0 * M_PI * t / n_sites * sin(2.0 * M_PI * (pick[k] - flux) / n_sites);
  }

  return s;
}

typedef struct {
  double re, im;
} cx;

static cx cmul(cx a, cx b) {
  return (cx){a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
static cx csub(cx a, cx b) { return (cx){a.re - b.re, a.im - b.im}; }
static double cabs2(cx a) { return a.re * a.re + a.im * a.im; }
static cx cdiv(cx a, cx b) {
  double d = cabs2(b);

  return (cx){(a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d};
}

// Solve A x = rhs for a dense complex system (augmented matrix a, n x (n + 1))
// by Gaussian elimination with partial pivoting
static int csolve(int n, cx a[AR_MAX][AR_MAX + 1], cx *x) {
  for (int c = 0; c < n; c++) {
    int p = c;

    for (int r = c + 1; r < n; r++) {
      if (cabs2(a[r][c]) > cabs2(a[p][c])) {
        p = r;
      }
    }

    if (cabs2(a[p][c]) < 1e-300) {
      return -1;
    }

    if (p != c) {
      for (int k = 0; k <= n; k++) {
        cx tmp = a[c][k];

        a[c][k] = a[p][k];
        a[p][k] = tmp;
      }
    }

    for (int r = c + 1; r < n; r++) {
      cx f = cdiv(a[r][c], a[c][c]);

      for (int k = c; k <= n; k++) {
        a[r][k] = csub(a[r][k], cmul(f, a[c][k]));
      }
    }
  }

  for (int r = n - 1; r >= 0; r--) {
    cx s = a[r][n];

    for (int k = r + 1; k < n; k++) {
      s = csub(s, cmul(a[r][k], x[k]));
    }

    x[r] = cdiv(s, a[r][r]);
  }

  return 0;
}

int ar_scatter(int n_sites, int a, int b, double flux, double energy, double t,
               ar_scatter_t *out) {
  if (n_sites < 4 || n_sites > AR_MAX || a < 0 || b < 0 || a >= n_sites ||
      b >= n_sites || a == b || !out || !isfinite(flux) || !isfinite(energy) ||
      !(t > 0.0) || !(fabs(energy) < 2.0 * t)) {
    return -1;
  }

  double gamma = sqrt(4.0 * t * t - energy * energy);
  cx sigma = {0.5 * energy, -0.5 * gamma}; // (E - i \Gamma) / 2
  double c = cos(2.0 * M_PI * flux), sn = sin(2.0 * M_PI * flux);
  static cx m[AR_MAX][AR_MAX + 1]; // augmented matrix (E - H - \Sigma | e_a)
  cx x[AR_MAX];

  for (int i = 0; i < n_sites; i++) {
    for (int j = 0; j <= n_sites; j++) {
      m[i][j] = (cx){0.0, 0.0};
    }

    m[i][i] = (cx){energy, 0.0};
  }

  for (int i = 0; i < n_sites; i++) {
    int j = (i + 1) % n_sites;
    cx h = (i == n_sites - 1) ? (cx){-t * c, -t * sn}
                              : (cx){-t, 0.0}; // flux phase on bond N-1 -> 0

    m[i][j] = csub(m[i][j], h);                 // H_ij = h
    m[j][i] = csub(m[j][i], (cx){h.re, -h.im}); // H_ji = conj(h)
  }

  m[a][a] = csub(m[a][a], sigma);
  m[b][b] = csub(m[b][b], sigma);
  m[a][n_sites] = (cx){1.0, 0.0};
  if (csolve(n_sites, m, x) != 0) {
    return -1;
  }

  // x = G e_a: G_ba = x[b], G_aa = x[a]
  out->transmission = gamma * gamma * cabs2(x[b]);
  cx saa = {1.0 + gamma * x[a].im, -gamma * x[a].re}; // 1 - i \Gamma G_aa

  out->reflection = cabs2(saa);

  return 0;
}
