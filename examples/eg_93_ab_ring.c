/*
 * Aharonov-Bohm ring
 *
 * NOTE: Thread magnetic flux through a tight-binding ring and every level
 * shifts although no field touches sites. Energies are periodic in  flux
 * quantum, so ground-state current is a periodic function of flux with sign set
 * by parity of electron number. Attach two leads at opposite sites and
 * transmission oscillates between a maximum at zero flux and exactly zero at
 * half a flux quantum: two arms interfere with a phase 2 \pi \Phi / \Phi_0
 */

#include "../physics/ab_ring.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  const int n = 16;

  printf(" > Persistent current of the N = %d ring (units of t / flux "
         "quantum), t = 1\n\n  %6s %10s %10s\n",
         n, "phi", "ne = 7", "ne = 8");
  for (int i = 1; i <= 9; i++) {
    double phi = 0.05 * i;

    printf("  %6.2f %10.5f %10.5f\n", phi, ar_current(n, 7, phi, 1.0),
           ar_current(n, 8, phi, 1.0));
  }

  printf(" > Transmission through leads at sites 0 and N/2 at E = 0.6\n\n  "
         "%6s %10s %10s\n",
         "phi", "T", "R");
  for (int i = 0; i <= 10; i++) {
    ar_scatter_t s;
    double phi = 0.05 * i;

    ar_scatter(n, 0, n / 2, phi, 0.6, 1.0, &s);
    printf("  %6.2f %10.5f %10.5f\n", phi, s.transmission, s.reflection);
  }

  return 0;
}
