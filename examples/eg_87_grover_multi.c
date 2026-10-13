/*
 * Grover search with several marked items
 *
 * NOTE: A database of N = 1024 items has M marked entries. Each Grover
 * iteration (oracle sign flip, inversion about mean) rotates state by 2 \theta,
 * \sin(\theta) = \sqrt(M / N), towards uniform superposition of marked items
 * The success probability is \sin^2((2k + 1) \theta): it peaks at k* ~ (\pi /
 * 4)
 * \sqrt(N / M) and then falls again if you keep iterating
 */

#include "../physics/grover_multi.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  const int n = 10;
  static double amp[1 << 10];
  static unsigned char marked[1 << 10];

  printf(" > Grover search, N = %d items\n\n", 1 << n);
  for (int m = 1; m <= 16; m *= 4) {
    for (int i = 0; i < (1 << n); i++) {
      marked[i] = (unsigned char)(i < m * 5 && i % 5 == 0);
    }

    grover_multi_init(amp, n);
    int kopt = grover_multi_optimal_iterations(n, m);

    printf("  M = %2d  \\theta = %.5f  best k = %d\n", m,
           grover_multi_theta(n, m), kopt);
    printf("  %4s  %12s  %12s\n", "k", "simulated", "\\sin^2");
    for (int k = 0; k <= 2 * kopt; k++) {
      if (k % ((kopt + 5) / 6) == 0 || k == kopt) {
        printf("  %4d  %12.6f  %12.6f%s\n", k,
               grover_multi_success(amp, n, marked),
               grover_multi_success_exact(n, m, k),
               k == kopt ? "  <- optimum" : "");
      }

      grover_multi_step(amp, n, marked);
    }

    printf("\n");
  }

  return 0;
}
