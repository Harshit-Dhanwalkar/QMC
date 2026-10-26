/*
 * Tunnelling in a double well
 *
 * NOTE: In a symmetric double well two lowest levels form a parity doublet
 * whose splitting \Delta shrinks exponentially with barrier. A particle placed
 * in left well, (\psi0 - \psi1) / \sqrt(2), tunnels to right and back with
 * period 2 \pi / \Delta: physics of ammonia inversion that drove first maser
 * Units \hbar = m = 1
 */

#include "../physics/double_well.h"
#include <math.h>
#include <stdio.h>

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
  enum { N = 1200 };
  static double v[N], psi[4 * N], re[N], im[N];
  const double hw = 5.0;
  const double a = 1.5;
  const double dx = dw_spacing(N, hw);
  double e[4];

  printf(" > Splitting of lowest doublet versus barrier height "
         "(well centres at +-%.1f)\n\n  %8s %14s %14s\n",
         a, "B", "E0", "\\Delta = E1 - E0");
  for (int b = 4; b <= 28; b += 4) {
    dw_potential(b, a, hw, N, v);
    dw_eigen(v, N, dx, 4, e, psi);

    printf("  %8d %14.6f %14.4e\n", b, e[0], e[1] - e[0]);
  }

  dw_potential(8.0, a, hw, N, v);
  dw_eigen(v, N, dx, 4, e, psi);

  const double c[4] = {M_SQRT1_2, -M_SQRT1_2, 0.0, 0.0};
  double d = e[1] - e[0];

  printf("\n > Packet prepared in left well at B = 8: "
         "\\Delta = %.5f, period 2 \\pi / \\Delta = %.1f\n\n"
         "  %10s %12s %10s\n",
         d, 2 * M_PI / d, "t / period", "P(left)", "<x>");
  for (int k = 0; k <= 8; k++) {
    double t = 2.0 * M_PI / d * k / 8.0;

    dw_evolve(psi, e, c, N, 2, t, re, im);
    printf("  %10.3f %12.5f %10.4f\n", k / 8.0, dw_prob_left(re, im, N, dx),
           dw_mean_x(re, im, N, dx, hw));
  }

  return 0;
}
