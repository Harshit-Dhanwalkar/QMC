/*
 * Quantum Zeno effect
 *
 * NOTE: A resonant drive flips a qubit in time T = \pi / \Omega
 * Measuring it N times along way suppresses flip (pulsed Zeno): survival -> 1
 * as N grows A continuous measurement, modelled as dephasing at rate \Gamma,
 * gives a decay rate that peaks at \Gamma = 2 \Omega and then falls as \Omega^2
 * / \Gamma: harder you look, slower it moves
 */

#include "../physics/zeno.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void) {
  const double omega = 1.0, T = M_PI / omega;

  printf(" > Pulsed measurements, T = \\pi / \\Omega\n\n  %6s  %12s\n", "N",
         "survival");
  for (int n = 1; n <= 256; n *= 2) {
    printf("  %6d  %12.6f\n", n, zeno_survival_pulsed(omega, T, n));
  }
  printf("\n > Continuous measurement: decay rate / Omega versus dephasing "
         "\\Gamma / \\Omega\n\n  %8s  %10s\n",
         "\\Gamma", "rate");

  for (int i = 0; i <= 12; i++) {
    double g = 0.5 * i;

    printf("  %8.2f  %10.5f\n", g, zeno_continuous_rate(omega, g));
  }

  printf("\n  Survival at t = %.3f with \\Gamma = 20 Omega: %.5f\n", T,
         zeno_survival_continuous(omega, 20.0, T));

  return 0;
}
