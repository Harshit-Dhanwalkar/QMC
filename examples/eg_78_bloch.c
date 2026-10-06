/*
 * A driven, damped qubit on the Bloch sphere
 *
 * NOTE: Units \hbar = 1, basis |g> = north pole (z = +1), |e> = south pole.
 * bloch_evolve() integrates the Lindblad equation (RK4) for
 *     H = (1/2)(\Omega \sigma_x + \Delta \sigma_z),
 *     L1 = \sqrt(\gamma1) |g><e|,
 *     L2 = \sqrt(\gamma_phi/2) \sigma_z
 * and every result is checked against a closed form:
 *
 *  1. Resonant Rabi oscillation: P_e(t) = \sin^2(\Omega t/2)
 *  2. \pi pulse: \Omega(t) = \pi flips |g> to |e>
 *  3. T1 relaxation of |e>: z(t) = 1 - 2 \exp^{-\gamma1 t}
 *  4. Dephasing of |+>: x(t) = \exp^{-\gamma_{\phi} t}
 *  5. Driven + damped steady state: P_e = (\Omega^2/4) /
 *     (\Delta^2 + \gamma1^2/4 + \Omega^2/2)
 */

#include "../physics/bloch.h"
#include "../physics/rabi.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
  printf(" > Bloch-sphere qubit dynamics\n\n");

  printf("  === Rabi oscillation (\\Omega = 1, resonant) ===\n");
  printf("  %6s %12s %12s\n", "t", "P_e (RK4)", "\\sin^2(t/2)");
  for (int k = 1; k <= 6; k++) {
    double v[3] = {0.0, 0.0, 1.0};
    double t = k * M_PI / 3.0;

    int n = (int)lround(t / 0.002);

    bloch_evolve(v, 1.0, 0.0, 0.0, 0.0, t / n, n); // dt = t/n hits t exactly
    printf("  %6.3f %12.6f %12.6f\n", t, 0.5 * (1.0 - v[2]),
           rabi_excited_probability(t, 1.0, 0.0));
  }

  printf("  === pi pulse (\\Omega = 2, t = \\pi/2) ===\n");
  double p[3] = {0.0, 0.0, 1.0};
  bloch_evolve(p, 2.0, 0.0, 0.0, 0.0, (M_PI / 2.0) / 800.0, 800);
  printf("  Bloch vector after the pulse: (%+.5f, %+.5f, %+.5f)\n", p[0], p[1],
         p[2]);

  printf("  === T1 relaxation of |e> (\\gamma1 = 0.5) ===\n");
  printf("  %6s %12s %12s\n", "t", "z (RK4)", "1-2\\exp^{-g1 t}");
  for (int k = 1; k <= 4; k++) {
    double v[3] = {0.0, 0.0, -1.0};
    double t = 1.0 * k;

    bloch_evolve(v, 0.0, 0.0, 0.5, 0.0, 0.002, (int)lround(t / 0.002));
    printf("  %6.2f %12.6f %12.6f\n", t, v[2], 1.0 - 2.0 * exp(-0.5 * t));
  }

  printf("  === Dephasing of |+> (\\gamma_{\\phi} = 0.4) ===\n");
  double d[3] = {1.0, 0.0, 0.0};
  bloch_evolve(d, 0.0, 0.0, 0.0, 0.4, 0.002, 1500); // t = 3
  printf("  x(3) = %.6f   e^{-gamma_phi t} = %.6f   z = %.1e\n", d[0],
         exp(-0.4 * 3.0), d[2]);

  printf("  === Steady state (\\Omega = 1.5, \\Delta = 0.5, \\gamma1 = 0.5) "
         "===\n");
  double s[3] = {0.0, 0.0, 1.0};
  bloch_evolve(s, 1.5, 0.5, 0.5, 0.0, 0.005, 10000); // t = 50
  double pe = 0.25 * 1.5 * 1.5 / (0.25 + 0.0625 + 0.5 * 2.25);
  printf("  P_e = %.6f   closed form = %.6f   |v| = %.4f (mixed)\n",
         0.5 * (1.0 - s[2]), pe, sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]));

  return 0;
}
