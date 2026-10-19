/*
 * Quantum corral: eigenmodes of a circular well
 *
 * NOTE: A particle in a unit disc has eigenfunctions J_m(k r) \cos(m \theta)
 * and J_m(k r) \sin(m \theta) with k zeros of Bessel function J_m, E = k^2
 * Every level with m >= 1 is doubly degenerate. This prints lowest levels with
 * their quantum numbers and compares number of modes below k with Weyl's law,
 * N(k) ~ k^2 / 4 - k / 2 (area term minus perimeter term)
 */

#include "../physics/circular_billiard.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  static cb_mode_t modes[512];
  int count = cb_modes(25.0, modes, 512);

  printf(" > %d modes with k <= 25\n\n  %4s  %3s  %3s  %5s  %10s\n", count, "#",
         "m", "n", "trig", "E = k^2");
  for (int i = 0; i < 16; i++) {
    printf("  %4d  %3d  %3d  %5s  %10.5f\n", i, modes[i].m, modes[i].n,
           modes[i].m == 0 ? "-" : (modes[i].sine ? "sin" : "cos"),
           modes[i].k * modes[i].k);
  }
  printf("\n  %6s  %8s  %10s\n", "k", "modes", "Weyl");
  for (int k = 5; k <= 25; k += 5) {
    int n = 0;

    for (int i = 0; i < count; i++) {
      n += modes[i].k <= k;
    }

    printf("  %6d  %8d  %10.1f\n", k, n, cb_weyl_count(k));
  }

  return 0;
}
