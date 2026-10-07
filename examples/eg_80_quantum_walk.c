/*
 * Discrete-time quantum walk versus classical random walk
 *
 * NOTE: Both walkers start at origin and take one step per tick
 * Classical width is \sqrt(t). With a Hadamard coin quantum walker
 * spreads as sqrt(1 - 1 / \sqrt(2)) t = 0.541 t, with probability piled up
 * near x = +/- t / \sqrt(2)
 */

#include "../physics/quantum_walk.h"
#include "vector.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { N = 2048 };

int main(void) {
  cvector_t *up = cvector_alloc(N);
  cvector_t *down = cvector_alloc(N);

  if (!up || !down) {
    return 1;
  }

  printf(" > Quantum walk (Hadamard coin, symmetric start)\n\n");
  printf("  %6s  %12s  %12s  %10s\n", "t", "quantum std", "classical", "ratio");
  
  qwalk_init(up, down, QWALK_COIN_SYM);

  int t = 0;
  for (int target = 10; target <= 640; target *= 2) {
    if (qwalk_step(up, down, M_PI / 4, target - t) != 0) {
      return 1;
    }

    t = target;
    double q = sqrt(qwalk_variance(up, down));
    double c = sqrt((double)t);

    printf("  %6d  %12.3f  %12.3f  %10.2f\n", t, q, c, q / c);
  }

  printf("\n  Asymptotic quantum width / t = %.4f\n",
         sqrt(qwalk_asymptotic_variance(M_PI / 4)));

  cvector_free(up);
  cvector_free(down);

  return 0;
}
