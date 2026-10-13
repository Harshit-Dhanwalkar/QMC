/*
 * CHSH Bell test on a noisy singlet
 *
 * NOTE: Two qubits share a Werner state, singlet mixed with white noise of
 * visibility V. At standard settings (Alice 0 and \pi/2, Bob pi/4 and 3 \pi/4)
 * quantum mechanics gives |S| = 2 \sqrt{2} V, above local-hidden-variable bound
 * 2 whenever V > 1 / \sqrt{2} = 0.7071. A concrete local model (a hidden angle,
 * each side outputs sign of projection) reaches exactly 2 at same settings
 */

#include "../physics/bell_chsh.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
  const double a = 0.0, a2 = M_PI / 2.0, b = M_PI / 4.0, b2 = 3.0 * M_PI / 4.0;

  printf(" > CHSH test, settings 0, 90, 45, 135 degrees\n\n");
  printf("  local-model S: %.4f   deterministic bound: %.0f\n\n",
         fabs(bell_local_chsh(a, a2, b, b2)), bell_deterministic_bound());
  printf("  %10s  %12s  %12s  %10s\n", "visibility", "|S| quantum",
         "grid maximum", "violates?");

  for (int i = 0; i <= 10; i++) {
    double v = 0.5 + 0.05 * i;
    double s = fabs(bell_chsh(v, a, a2, b, b2));

    printf("  %10.2f  %12.5f  %12.5f  %10s\n", v, s, bell_max_chsh(v, 8),
           s > 2.0 ? "yes" : "no");
  }

  double sums[4];
  double est = bell_sample(1.0, a, a2, b, b2, 100000, 7u, sums);

  printf(
      "\n  100000 simulated pairs per setting, V = 1: S = %.4f (exact %.4f)\n",
      fabs(est), 2.0 * M_SQRT2);

  return 0;
}
