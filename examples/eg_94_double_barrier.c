/*
 * Resonant tunnelling through a double barrier
 *
 * NOTE: Two thin barriers around a well transmit almost nothing at most
 * energies, but at quasi-bound levels of well transmission jumps to exactly 1
 * and particle dwells there for a long time. This is Fabry-Perot effect of
 * matter waves, principle of resonant-tunnelling diode
 * Units \hbar^2 / 2m = 1, so E = k^2
 */

#include "../physics/double_barrier.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
  const double v[3] = {30.0, 0.0, 30.0};
  const double d[3] = {0.5, 3.0, 0.5};

  printf(" > Barriers V = 30, thickness 0.5, well width 3 (infinite-well "
         "levels %.3f, %.3f, %.3f)\n\n  %8s %12s %12s\n",
         pow(M_PI / 3, 2), pow(2 * M_PI / 3, 2), pow(3 * M_PI / 3, 2), "E", "T",
         "dwell");
  for (int i = 1; i <= 24; i++) {
    double e = 0.5 * i;

    printf("  %8.2f %12.6f %12.4f\n", e, db_transmission(3, v, d, e),
           db_dwell_time(3, v, d, e));
  }

  printf("\n > Single barrier of same height and thickness at same energies: "
         "T(E=2.2) = %.6f\n",
         db_transmission(1, v, d, 2.2));

  return 0;
}
