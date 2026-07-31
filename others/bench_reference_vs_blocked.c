/*
Benchmark: reference (unblocked) vs column-blocked tridiagonal QL eigenvector
computation.

both implementations are mathematically identical

The "blocked" variant tests a fix for that record every Givens rotation during
(cheap, cache-friendly) eigenvalue-only pass, then replay them against
eigenvector matrix one column-block at a time, so working set per pass is
block_size*n instead of n*n.

USAGE:
    gcc -O2 -o bench_reference_vs_blocked bench_reference_vs_blocked.c -lm
    ./bench_reference_vs_blocked <N> <mode>
        mode = 0        -> reference (unblocked) implementation
        mode = <block>  -> blocked implementation with that column-block
                            size (e.g. 32, 64, 128, 256)

Suggested runs to reproduce what was measured in sandbox, and to get real
diagnostics locally: # Wall-clock comparison at the sizes where anomaly was
worst for mode in 0 64 128 256; do ./bench_reference_vs_blocked 1600 $mode; done
  for mode in 0 64 128 256; do ./bench_reference_vs_blocked 3200 $mode; done

 # Cachegrind
 valgrind --tool=cachegrind --cache-sim=yes ./bench_reference_vs_blocked 3200 0
 valgrind --tool=cachegrind --cache-sim=yes ./bench_reference_vs_blocked 3200
128 # then: cg_annotate cachegrind.out.<pid>
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double pythag(double a, double b) {
  double absa = fabs(a), absb = fabs(b);
  if (absa > absb) {
    double r = absb / absa;
    return absa * sqrt(1.0 + r * r);
  } else if (absb == 0.0)
    return 0.0;
  else {
    double r = absa / absb;
    return absb * sqrt(1.0 + r * r);
  }
}

typedef struct {
  int i;
  double s, c;
} rot_t;

// Reference (original, unblocked) implementation for correctness comparison
static void reference_eigh(const double *diag, const double *offdiag, int n,
                           double *d, double *z) {
  double *e = malloc(n * sizeof *e);
  memcpy(d, diag, n * sizeof(double));
  for (int i = 0; i < n - 1; i++) {
    e[i] = offdiag[i];
  }

  e[n - 1] = 0.0;
  for (int ii = 0; ii < n; ii++) {
    for (int jj = 0; jj < n; jj++) {
      z[ii * n + jj] = (ii == jj) ? 1.0 : 0.0;
    }
  }

  for (int l = 0; l < n; l++) {
    int iter = 0, m;
    do {
      for (m = l; m < n - 1; m++) {
        double dd = fabs(d[m]) + fabs(d[m + 1]);
        if (fabs(e[m]) <= dd * 1e-15) {
          break;
        }
      }

      if (m != l) {
        if (++iter > 100) {
          break;
        }

        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        double r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? fabs(r) : -fabs(r)));
        double s = 1.0, c = 1.0, p = 0.0;
        int i;
        for (i = m - 1; i >= l; i--) {
          double f = s * e[i];
          double b = c * e[i];
          r = pythag(f, g);
          e[i + 1] = r;

          if (r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;

            break;
          }

          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;

          for (int k = 0; k < n; k++) {
            double f2 = z[(i + 1) * n + k];
            z[(i + 1) * n + k] = s * z[i * n + k] + c * f2;
            z[i * n + k] = c * z[i * n + k] - s * f2;
          }
        }
        if (r == 0.0 && i >= l) {
          continue;
        }

        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }
  free(e);
}

// Blocked implementation
static rot_t *g_rots;
static long g_nrots, g_cap;

static void record_rot(int i, double s, double c) {
  if (g_nrots >= g_cap) {
    g_cap = g_cap ? g_cap * 2 : 1024;
    g_rots = realloc(g_rots, g_cap * sizeof(rot_t));
  }

  g_rots[g_nrots].i = i;
  g_rots[g_nrots].s = s;
  g_rots[g_nrots].c = c;
  g_nrots++;
}

static void blocked_eigh(const double *diag, const double *offdiag, int n,
                         double *d, double *z, int block_size) {
  double *e = malloc(n * sizeof *e);
  memcpy(d, diag, n * sizeof(double));

  for (int i = 0; i < n - 1; i++) {
    e[i] = offdiag[i];
  }
  e[n - 1] = 0.0;

  g_nrots = 0;
  g_cap = 0;
  g_rots = NULL;

  // Phase 1: eigenvalues only, recording rotations
  for (int l = 0; l < n; l++) {
    int iter = 0, m;
    do {
      for (m = l; m < n - 1; m++) {
        double dd = fabs(d[m]) + fabs(d[m + 1]);

        if (fabs(e[m]) <= dd * 1e-15) {
          break;
        }
      }

      if (m != l) {
        if (++iter > 100) {
          break;
        }

        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        double r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? fabs(r) : -fabs(r)));
        double s = 1.0, c = 1.0, p = 0.0;
        int i;

        for (i = m - 1; i >= l; i--) {
          double f = s * e[i];
          double b = c * e[i];
          r = pythag(f, g);
          e[i + 1] = r;

          if (r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;

            break;
          }

          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          record_rot(i, s, c);
        }

        if (r == 0.0 && i >= l) {
          continue;
        }
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }
  free(e);

  // Phase 2: initialize z to identity, replay rotations column-block by
  // column-block
  for (int ii = 0; ii < n; ii++) {
    for (int jj = 0; jj < n; jj++) {
      z[ii * n + jj] = (ii == jj) ? 1.0 : 0.0;
    }
  }

  for (int k0 = 0; k0 < n; k0 += block_size) {
    int k1 = k0 + block_size;
    if (k1 > n) {
      k1 = n;
    }

    for (long r = 0; r < g_nrots; r++) {
      int i = g_rots[r].i;
      double s = g_rots[r].s, c = g_rots[r].c;
      for (int k = k0; k < k1; k++) {
        double f2 = z[(i + 1) * n + k];
        z[(i + 1) * n + k] = s * z[i * n + k] + c * f2;
        z[i * n + k] = c * z[i * n + k] - s * f2;
      }
    }
  }
  free(g_rots);
}

int main(int argc, char **argv) {
  int n = atoi(argv[1]);
  int mode = atoi(argv[2]); // 0=reference, else=block size
  double *diag = malloc(n * sizeof(double)),
         *offdiag = malloc((n - 1) * sizeof(double));

  for (int i = 0; i < n; i++) {
    diag[i] = 2.0 + 0.001 * i;
  }
  for (int i = 0; i < n - 1; i++) {
    offdiag[i] = -1.0;
  }

  double *d = malloc(n * sizeof(double)),
         *z = malloc((size_t)n * n * sizeof(double));

  clock_t t0 = clock();
  if (mode == 0) {
    reference_eigh(diag, offdiag, n, d, z);
  } else {
    blocked_eigh(diag, offdiag, n, d, z, mode);
  }

  clock_t t1 = clock();
  printf("N=%5d mode=%-6d time=%8.4fs\n", n, mode,
         (double)(t1 - t0) / CLOCKS_PER_SEC);
  fflush(stdout);

  return 0;
}
