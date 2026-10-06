/*
 * Hydrogen orbitals by Monte Carlo sampling
 *
 * NOTE: Draws points from |\psi_{nlm}|^2 (radial inverse-CDF of analytic
 * R_{nl}, rejection-sampled spherical harmonics) and compares sample averages
 * with exact hydrogen results
 *     <r>   = (3 n^2 - l(l+1)) / 2           (in Bohr radii)
 *     <r^2> = n^2 (5 n^2 + 1 - 3 l(l+1)) / 2
 * shows how real p_z orbital splits its density between two lobes (phase 0
 * above xy plane, \pi below)
 */

#include "../physics/orbital_sample.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Point i of a packed (x, y, z, phase) array; size_t arithmetic avoids int
// overflow
static const double *pt(const double *p, int i) { return p + 4 * (size_t)i; }

enum { COUNT = 100000 };

int main(void) {
  printf(" > Hydrogen orbital sampling (%d points each)\n\n", (int)COUNT);

  double *pts = malloc(4 * (size_t)COUNT * sizeof *pts);
  if (!pts) {
    return 1;
  }

  const int nl[][2] = {{1, 0}, {2, 0}, {2, 1}, {3, 2}, {4, 3}};

  printf("  %-6s %10s %10s %12s %12s\n", "orbital", "<r> MC", "<r> exact",
         "<r^2> MC", "<r^2> exact");
  for (int k = 0; k < 5; k++) {
    int n = nl[k][0];
    int l = nl[k][1];

    if (hydrogen_orbital_sample(n, l, 0, 0, COUNT, 12345 + k, pts) != 0) {
      free(pts);

      return 1;
    }

    double s1 = 0.0;
    double s2 = 0.0;
    for (int i = 0; i < COUNT; i++) {
      double r2 = pt(pts, i)[0] * pt(pts, i)[0] +
                  pt(pts, i)[1] * pt(pts, i)[1] + pt(pts, i)[2] * pt(pts, i)[2];

      s1 += sqrt(r2);
      s2 += r2;
    }

    printf("  %d%c     %10.4f %10.4f %12.3f %12.3f\n", n, "spdfg"[l],
           s1 / COUNT, (3.0 * n * n - l * (l + 1.0)) / 2.0, s2 / COUNT,
           n * n * (5.0 * n * n + 1.0 - 3.0 * l * (l + 1.0)) / 2.0);
  }

  // Real 2p_z: two lobes of opposite sign
  hydrogen_orbital_sample(2, 1, 0, 1, COUNT, 99, pts);
  int upper = 0;
  int upper_positive = 0;
  for (int i = 0; i < COUNT; i++) {
    if (pt(pts, i)[2] > 0.0) {
      upper++;
      upper_positive += (pt(pts, i)[3] < 1.0);
    }
  }

  printf("\n  2p_z: %.1f%% of dots are above xy plane, %.1f%% of those "
         "carry \\psi > 0\n",
         100.0 * upper / COUNT, 100.0 * upper_positive / (upper ? upper : 1));

  free(pts);

  return 0;
}
