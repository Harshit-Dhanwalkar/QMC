/*
 * Tunnelling in a symmetric double well (ammonia-inversion physics)
 *
 * Eigenstates by Sturm bisection + parity-projected inverse iteration; time
 * evolution by projection onto the low-lying eigenbasis
 */

#include "double_well.h"
#include <math.h>
#include <stddef.h>

double dw_spacing(int n, double half_width) {
  return 2.0 * half_width / (n + 1);
}

int dw_potential(double barrier, double a, double half_width, int n,
                 double *v) {
  if (!v || n < 3 || n > DW_MAX_GRID || !(barrier >= 0.0) || !(a > 0.0) ||
      !(half_width > 0.0) || !isfinite(barrier) || !isfinite(a) ||
      !isfinite(half_width)) {
    return -1;
  }

  double dx = dw_spacing(n, half_width);
  for (int i = 0; i < n; i++) {
    double x = -half_width + (i + 1) * dx;
    double u = 1.0 - x * x / (a * a);

    v[i] = barrier * u * u;
  }

  return 0;
}

// Number of eigenvalues of tridiagonal H below e (Sturm sequence)
static int count_below(const double *d, double off, int n, double e) {
  int cnt = 0;
  double q = d[0] - e;

  if (q < 0) {
    cnt++;
  }
  for (int i = 1; i < n; i++) {
    if (q == 0.0) {
      q = 1e-300;
    }

    q = d[i] - e - off * off / q;
    if (q < 0) {
      cnt++;
    }
  }

  return cnt;
}

// Solve (H - shift) x = b for tridiagonal H (constant off-diagonal) by
// Thomas algorithm; work holds 2 n doubles
static void tri_solve(const double *d, double off, int n, double shift,
                      const double *b, double *x, double *work) {
  double *cp = work, *dp = work + n;
  double denom = d[0] - shift;

  if (fabs(denom) < 1e-300) {
    denom = 1e-300;
  }

  cp[0] = off / denom;
  dp[0] = b[0] / denom;
  for (int i = 1; i < n; i++) {
    denom = (d[i] - shift) - off * cp[i - 1];
    if (fabs(denom) < 1e-300) {
      denom = 1e-300;
    }

    cp[i] = off / denom;
    dp[i] = (b[i] - off * dp[i - 1]) / denom;
  }

  x[n - 1] = dp[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    x[i] = dp[i] - cp[i] * x[i + 1];
  }
}

int dw_eigen(const double *v, int n, double dx, int n_states, double *e,
             double *psi) {
  if (!v || !e || !psi || n < 3 || n > DW_MAX_GRID || n_states < 1 ||
      n_states > DW_MAX_STATES || n_states > n || !(dx > 0.0)) {
    return -1;
  }

  double scale = 0.0;

  for (int i = 0; i < n; i++) {
    scale = fmax(scale, fabs(v[i]));
  }

  for (int i = 0; i < n / 2; i++) {
    if (fabs(v[i] - v[n - 1 - i]) > 1e-9 * (scale + 1.0)) {
      return -1;
    }
  }

  static double d[DW_MAX_GRID], work[2 * DW_MAX_GRID], rhs[DW_MAX_GRID],
      x[DW_MAX_GRID];
  double off = -0.5 / (dx * dx);
  double lo = 1e300;
  double hi = -1e300;

  for (int i = 0; i < n; i++) {
    d[i] = 1.0 / (dx * dx) + v[i];
    lo = fmin(lo, d[i] - fabs(off) * 2.0);
    hi = fmax(hi, d[i] + fabs(off) * 2.0);
  }

  for (int k = 0; k < n_states; k++) {
    // k-th eigenvalue (0-based): smallest e with count_below(e) >= k + 1
    double a = lo;
    double b = hi;

    for (int it = 0; it < 200 && b - a > 1e-15 * (fabs(a) + fabs(b) + 1e-12);
         it++) {
      double m = 0.5 * (a + b);

      if (count_below(d, off, n, m) >= k + 1) {
        b = m;
      } else {
        a = m;
      }
    }

    e[k] = 0.5 * (a + b);

    // parity-projected inverse iteration
    double *vec = &psi[(size_t)k * n];
    int parity = (k % 2 == 0) ? 1 : -1;
    unsigned seed = 12345u + 7919u * (unsigned)k;

    for (int i = 0; i < n; i++) {
      seed = seed * 1664525u + 1013904223u;
      vec[i] = (double)(seed >> 8) / 16777216.0 - 0.5;
    }

    for (int it = 0; it < 8; it++) {
      // project onto parity sector
      for (int i = 0; i < n / 2; i++) {
        double s = 0.5 * (vec[i] + parity * vec[n - 1 - i]);

        vec[i] = s;
        vec[n - 1 - i] = parity * s;
      }

      if (n % 2) {
        vec[n / 2] = parity > 0 ? vec[n / 2] : 0.0;
      }

      for (int i = 0; i < n; i++) {
        rhs[i] = vec[i];
      }

      tri_solve(d, off, n, e[k] + 1e-13 * (fabs(e[k]) + 1.0), rhs, x, work);
      double nrm = 0.0;

      for (int i = 0; i < n; i++) {
        nrm += x[i] * x[i];
      }

      nrm = sqrt(nrm);
      for (int i = 0; i < n; i++) {
        vec[i] = x[i] / nrm;
      }

      // Gram-Schmidt against lower states of same parity
      for (int j = 0; j < k; j++) {
        const double *o = &psi[(size_t)j * n];
        double dot = 0.0;

        for (int i = 0; i < n; i++) {
          dot += vec[i] * o[i];
        }
        for (int i = 0; i < n; i++) {
          vec[i] -= dot * o[i];
        }
      }
      nrm = 0.0;
      for (int i = 0; i < n; i++) {
        nrm += vec[i] * vec[i];
      }
      nrm = sqrt(nrm);
      for (int i = 0; i < n; i++) {
        vec[i] /= nrm;
      }
    }

    // normalise on grid and fix sign
    double sum = 0.0;
    double best = 0.0;
    int bi = 0;

    for (int i = 0; i < n; i++) {
      sum += vec[i] * vec[i] * dx;
    }

    sum = sqrt(sum);
    // even: largest |psi| anywhere; odd: largest |\psi| on right half
    for (int i = 0; i < n; i++) {
      vec[i] /= sum;
      if ((parity > 0 || i > n / 2) && fabs(vec[i]) > best) {
        best = fabs(vec[i]);
        bi = i;
      }
    }

    if (vec[bi] < 0.0) {
      for (int i = 0; i < n; i++) {
        vec[i] = -vec[i];
      }
    }
  }

  return 0;
}

int dw_project(const double *psi, int n, double dx, int n_states,
               double half_width, double center, double sigma, double *c,
               double *captured) {
  if (!psi || !c || !captured || n < 3 || n > DW_MAX_GRID || n_states < 1 ||
      n_states > DW_MAX_STATES || !(sigma > 0.0) || !(dx > 0.0) ||
      !isfinite(center)) {
    return -1;
  }

  static double g[DW_MAX_GRID];
  double norm = 0.0;

  for (int i = 0; i < n; i++) {
    double x = -half_width + (i + 1) * dx;

    g[i] = exp(-(x - center) * (x - center) / (4.0 * sigma * sigma));
    norm += g[i] * g[i] * dx;
  }

  norm = sqrt(norm);
  *captured = 0.0;
  for (int k = 0; k < n_states; k++) {
    double s = 0.0;

    for (int i = 0; i < n; i++) {
      s += psi[(size_t)k * n + i] * g[i] / norm * dx;
    }

    c[k] = s;
    *captured += s * s;
  }

  return 0;
}

int dw_evolve(const double *psi, const double *e, const double *c, int n,
              int n_states, double t, double *re, double *im) {
  if (!psi || !e || !c || !re || !im || n < 1 || n_states < 1 || !isfinite(t)) {
    return -1;
  }

  for (int i = 0; i < n; i++) {
    re[i] = im[i] = 0.0;
  }

  for (int k = 0; k < n_states; k++) {
    double cr = c[k] * cos(e[k] * t);
    double ci = -c[k] * sin(e[k] * t);

    for (int i = 0; i < n; i++) {
      re[i] += cr * psi[(size_t)k * n + i];
      im[i] += ci * psi[(size_t)k * n + i];
    }
  }

  return 0;
}

double dw_prob_left(const double *re, const double *im, int n, double dx) {
  double left = 0.0;
  double tot = 0.0;

  for (int i = 0; i < n; i++) {
    double p = re[i] * re[i] + im[i] * im[i];

    tot += p;
    if (i < n / 2) {
      left += p;
    }
  }

  (void)dx;

  return tot > 0.0 ? left / tot : NAN;
}

double dw_mean_x(const double *re, const double *im, int n, double dx,
                 double half_width) {
  double s = 0.0;
  double tot = 0.0;

  for (int i = 0; i < n; i++) {
    double p = re[i] * re[i] + im[i] * im[i];

    tot += p;
    s += p * (-half_width + (i + 1) * dx);
  }

  return tot > 0.0 ? s / tot : NAN;
}
