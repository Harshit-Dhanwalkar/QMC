/*
Instrumentation: counts outer QL sweeps and inner rotation ("i-loop") steps for
tridiagonal QL algorithm (eigenvalues only), directly on same clustered-diagonal
test matrix used in wall-clock benchmarks, across increasing N.

USAGE:
    gcc -O2 -o instrument_iteration_counts instrument_iteration_counts.c -lm
    for N in 400 800 1600 3200; do ./instrument_iteration_counts $N; done
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double pythag(double a, double b) {
  double absa = fabs(a), absb = fabs(b);
  if (absa > absb) {
    double r = absb / absa;
    return absa * sqrt(1.0 + r * r);
  } else if (absb == 0.0) {
    return 0.0;
  } else {
    double r = absa / absb;
    return absb * sqrt(1.0 + r * r);
  }
}

int main(int argc, char **argv) {
  int n = atoi(argv[1]);
  double *diag = malloc(n * sizeof(double)),
         *offdiag = malloc((n - 1) * sizeof(double));

  // exact same matrix as wall-clock benchmark: clustered diag
  for (int i = 0; i < n; i++) {
    diag[i] = 2.0 + 0.001 * i;
  }
  for (int i = 0; i < n - 1; i++) {
    offdiag[i] = -1.0;
  }

  double *d = malloc(n * sizeof(double)), *e = malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) {
    d[i] = diag[i];
  }
  for (int i = 0; i < n - 1; i++) {
    e[i] = offdiag[i];
  }
  e[n - 1] = 0.0;

  long total_outer = 0, total_i_steps = 0;
  int max_iter_seen = 0;

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
        total_outer++;
        if (++iter > 100) {
          break;
        }
        if (iter > max_iter_seen) {
          max_iter_seen = iter;
        }

        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        double r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? fabs(r) : -fabs(r)));
        double s = 1.0, c = 1.0, p = 0.0;
        int i;
        for (i = m - 1; i >= l; i--) {
          total_i_steps++;
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

  printf("N=%6d  total_outer=%10ld  total_i_steps=%12ld  i_steps/N^2=%.4f  "
         "i_steps/N^2.5=%.6f  i_steps/N^3=%.8f  max_iter_for_one_l=%d\n",
         n, total_outer, total_i_steps, (double)total_i_steps / ((double)n * n),
         (double)total_i_steps / pow(n, 2.5),
         (double)total_i_steps / ((double)n * n * n), max_iter_seen);

  return 0;
}
