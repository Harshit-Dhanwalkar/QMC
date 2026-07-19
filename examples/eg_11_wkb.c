/*
 * WKB Approximation
 *
 * Bohr-Sommerfeld quantization for harmonic oscillator (compared to exact E_n =
 * \hbar * \omega * (n+1/2)), and WKB tunneling through rectangular barrier
 * (compared to exact closed-form transmission coeff.).
 */

#include "../physics/wkb.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  printf(" > WKB Approximation\n\n");

  double hbar = 1.0, omega = 1.0;

  printf("   Bohr-Sommerfeld quantization, harmonic oscillator:\n");
  printf("   n    E_WKB     E_exact    error\n");
  printf("   --   -------   -------   --------\n");
  for (int n = 0; n < 5; n++) {
    double E_wkb = wkb_energy_harmonic(n, hbar, omega);
    double E_exact = hbar * omega * (n + 0.5);
    printf("   %2d   %6.4f   %6.4f   %7.2e\n", n, E_wkb, E_exact,
           fabs(E_wkb - E_exact));
  }
  printf("\n");

  printf("   WKB tunneling through a rectangular barrier:\n");
  printf("   E     V0    width   T_WKB (exact closed form)\n");
  printf("   ---   ---   -----   -------------------------\n");
  double hbar_sq_2m = 0.5;
  double cases[3][3] = {{0.3, 1.0, 1.0}, {0.5, 1.0, 2.0}, {0.2, 2.0, 1.5}};
  for (int i = 0; i < 3; i++) {
    double E = cases[i][0], V0 = cases[i][1], width = cases[i][2];
    double T = wkb_tunneling_rectangular(E, V0, width, hbar_sq_2m);
    printf("   %.1f   %.1f   %.1f     %.6e\n", E, V0, width, T);
  }

  return 0;
}
