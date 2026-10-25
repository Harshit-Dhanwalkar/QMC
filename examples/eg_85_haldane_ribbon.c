/*
 * Haldane model: chiral edge states on a ribbon
 *
 * NOTE: Graphene with complex next-nearest-neighbour hopping t2 \exp(+-i \phi)
 * and a sublattice mass M is a Chern insulator :
 *  (C = +-1) for |M| < 3 \sqrt(3) t2 |\sin(\phi)| and trivial otherwise
 * A 20-cell ribbon of Chern phase has one state per edge inside bulk gap
 * Prints states closest to E = 0 at a few wave numbers with their position
 * across ribbon (-1 first edge, +1 far edge), then same for a trivial mass
 */

#include "../physics/tight_binding.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { WIDTH = 20 };

static void scan(double mass) {
  double t2 = 0.1;
  double phi = M_PI / 2;
  double e[2 * WIDTH];
  double mc[2 * WIDTH];

  tb_model_t *m = tb_model_haldane(1.0, t2, phi, mass);
  if (!m) {
    return;
  }

  printf("\n  M = %.2f  (Chern number %d, bulk gap at K %.3f)\n", mass,
         tb_haldane_chern(1.0, t2, phi, mass),
         tb_haldane_valley_gap(t2, phi, mass, 1));
  printf("  %6s  %24s\n", "k1", "two states nearest E = 0");

  for (int i = 0; i <= 6; i++) {
    double k = M_PI * (0.5 + i / 6.0);

    tb_ribbon_bands(m, k, WIDTH, e, mc);

    printf("  %6.3f  E = %+7.4f (pos %+5.2f)   E = %+7.4f (pos %+5.2f)\n", k,
           e[WIDTH - 1], mc[WIDTH - 1], e[WIDTH], mc[WIDTH]);
  }

  tb_model_free(m);
}

int main(void) {
  printf(" > Haldane ribbon, %d cells, t1 = 1, t2 = 0.1, \\phi = \\pi/2\n",
         WIDTH);

  scan(0.0);
  scan(1.0);

  return 0;
}
