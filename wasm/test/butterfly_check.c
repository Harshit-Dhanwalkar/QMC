/* Native checks for wasm/qmc_butterfly.c, plus a digest that
 * wasm/test/butterfly_check.mjs must reproduce from the WASM build.
 *
 *   1. band edges from 4 special k-points == dense k-grid (all coprime p/q,
 *      q<=14)
 *   2. TKNN Hall numbers == library Fukui-Hatsugai-Suzuki Chern numbers
 *      (all coprime p/q, q<=9, every gapped r); middle gap of even q is closed
 *   3. digest line over all coprime p/q with q<=24 for native-vs-WASM diff
 */
#include "../../physics/tight_binding.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double *qmc_bf_buffer(void);
int qmc_bf_edges(int p, int q);
int qmc_bf_tknn(int p, int q, int r);
double qmc_bf_chern(int p, int q, int r, int n_k);

static int gcd(int a, int b) { return b ? gcd(b, a % b) : a; }

int main(void) {
  int fails = 0;

  // Band edges vs dense grid (k2 grid contains \pi/q exactly)
  double worst = 0.0;
  for (int q = 2; q <= 14; q++) {
    for (int p = 1; p < q; p++) {
      if (gcd(p, q) != 1) {
        continue;
      }

      tb_model_t *m = tb_model_hofstadter(p, q, 1.0);
      double lo[64], hi[64], e[64];
      for (int b = 0; b < q; b++) {
        lo[b] = 1e9;
        hi[b] = -1e9;
      }

      int n1 = 24, n2 = 8 * q;
      for (int a = 0; a < n1; a++)
        for (int c = 0; c < n2; c++) {
          tb_bands(m, 2 * M_PI * a / n1, 2 * M_PI * c / n2, e);
          for (int b = 0; b < q; b++) {
            if (e[b] < lo[b]) {
              lo[b] = e[b];
            }
            if (e[b] > hi[b]) {
              hi[b] = e[b];
            }
          }
        }

      tb_model_free(m);

      if (qmc_bf_edges(p, q) != 0) {
        printf("edges failed %d/%d\n", p, q);
        fails++;

        continue;
      }

      const double *g = qmc_bf_buffer();
      for (int b = 0; b < q; b++) {
        double d = fmax(fabs(lo[b] - g[2 * b]), fabs(hi[b] - g[2 * b + 1]));
        if (d > worst) {
          worst = d;
        }
      }
    }
  }

  printf("edges: worst |special-points - dense grid| = %.2e\n", worst);
  if (worst > 1e-9) {
    printf("FAIL: band edges\n");
    fails++;
  }

  // TKNN vs library Chern numbers
  int checked = 0, closed = 0;
  for (int q = 3; q <= 9; q++) {
    for (int p = 1; p < q; p++) {
      if (gcd(p, q) != 1) {
        continue;
      }

      for (int r = 1; r < q; r++) {
        int t = qmc_bf_tknn(p, q, r);
        double c = qmc_bf_chern(p, q, r, 16);
        if (isnan(c)) {
          closed++;
          if (t != INT_MIN) {
            printf("FAIL: gap closed but tknn=%d (p/q=%d/%d r=%d)\n", t, p, q,
                   r);
            fails++;
          }

          continue;
        }

        checked++;
        if (t == INT_MIN || fabs(c - t) > 0.05) {
          printf("FAIL: p/q=%d/%d r=%d tknn=%d numeric=%.3f\n", p, q, r, t, c);
          fails++;
        }
      }
    }
  }

  printf("chern: %d gaps agree with TKNN, %d closed gaps (even q middle)\n",
         checked, closed);

  // Digest
  double sum = 0.0;
  long tsum = 0;
  for (int q = 2; q <= 24; q++)
    for (int p = 1; p < q; p++) {
      if (gcd(p, q) != 1) {
        continue;
      }

      qmc_bf_edges(p, q);
      const double *g = qmc_bf_buffer();
      for (int b = 0; b < q; b++) {
        sum +=
            (g[2 * b + 1] - g[2 * b]) * (b + 1) + g[2 * b] * g[2 * b] * (q - b);
      }

      for (int r = 1; r < q; r++) {
        int t = qmc_bf_tknn(p, q, r);
        if (t != INT_MIN) {
          tsum += labs(t) * r;
        }
      }
    }

  printf("digest edges=%.9f tknn=%ld\n", sum, tsum);

  if (fails) {
    printf("%d FAILURES\n", fails);
    return 1;
  }
  printf("OK\n");
  return 0;
}
