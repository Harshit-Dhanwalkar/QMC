/*
 * Sequential Stern-Gerlach analyzers
 *
 * NOTE: A beam of spin-1/2 atoms is sent through analyzers with different axes.
 * Blocking one output is a projective measurement: z+ -> x+ -> z+ passes a
 * quarter of beam although z+ -> z+ passes all of it and z+ -> z- none
 * Passing both outputs through a recombining magnet restores state only to
 * extent that two packets overlap
 */

#include "../physics/stern_gerlach.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void show(const char *name, sg_vec_t r0, const sg_analyzer_t *a, int n) {
  double inten[8];
  sg_vec_t r[8];

  sg_chain(r0, 1.0, a, n, inten, r);
  printf("  %-44s", name);
  for (int k = 0; k <= n; k++) {
    printf(" %7.4f", inten[k]);
  }
  printf("\n");
}

int main(void) {
  sg_vec_t zup = {0, 0, 1};
  sg_analyzer_t zz[2] = {{0, 0, SG_PASS_UP, 1}, {0, 0, SG_PASS_DOWN, 1}};
  sg_analyzer_t zxz[3] = {{0, 0, SG_PASS_UP, 1},
                          {M_PI / 2, 0, SG_PASS_UP, 1},
                          {0, 0, SG_PASS_UP, 1}};
  sg_analyzer_t hd[2] = {{M_PI / 2, 0, SG_PASS_BOTH, 1.0},
                         {0, 0, SG_PASS_DOWN, 1}};

  printf(" > Intensity after each analyzer (source = 1)\n\n");
  show("z+ -> z-", zup, zz, 2);
  show("z+ -> x+ -> z+", zup, zxz, 3);
  show("z+ -> [x box, perfect] -> z-", zup, hd, 2);
  hd[0].coherence = 0.0;
  show("z+ -> [x box, no recombination] -> z-", zup, hd, 2);

  printf("\n > Rotating analyzer by pi / N each step (z+ in, z- out)\n\n  "
         "%4s %12s\n",
         "N", "transmitted");
  for (int n = 1; n <= 64; n *= 2) {
    sg_analyzer_t st[64];
    double inten[65];
    sg_vec_t r[65];

    for (int k = 0; k < n; k++) {
      st[k].theta = (k + 1) * M_PI / n;
      st[k].phi = 0.0;
      st[k].mode = SG_PASS_UP;
      st[k].coherence = 1.0;
    }
    sg_chain(zup, 1.0, st, n, inten, r);
    printf("  %4d %12.5f\n", n, inten[n]);
  }

  printf("\n > Packet overlap \\exp(-d^2 / 8 \\sigma^2) as magnet separates "
         "beams\n\n  %8s %10s\n",
         "d/\\sigma", "overlap");
  for (int d = 0; d <= 6; d++) {
    printf("  %8d %10.5f\n", d, sg_coherence(d, 1.0));
  }

  return 0;
}
