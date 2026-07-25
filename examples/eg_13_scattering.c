/*
 * Scattering Theory (Phase Shifts and Born Approximation)
 *
 * s-wave phase shift for a hard-sphere-like barrier, and Born
 * approximation differential cross section for a screened Coulomb
 * (Yukawa) potential.
 */

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/utils.h"
#include "../physics/potentials.h"
#include "../physics/scattering.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Scattering Theory: Phase Shifts and Born Approximation\n\n");

  double a = 2.0, V0 = 1000.0;
  double barrier_params[2] = {a, V0};
  double hbar_sq_2m = 0.5;

  printf("   s-wave phase shift, barrier radius a=%.1f (exact hard-sphere: "
         "\\delta_0=-ka)\n",
         a);
  printf("   k      \\delta_0 (num)   -ka (exact)\n");
  printf("   ----   -------------   -----------\n");
  for (double k = 0.2; k <= 1.0; k += 0.2) {
    double delta = phase_shift(0, k, V_barrier, barrier_params, 0.01, 20.0,
                               2000, hbar_sq_2m);
    printf("   %4.2f   %13.4f   %11.4f\n", k, delta, -k * a);
  }
  printf("\n");

  double g = 1.0, mu = 1.0;
  double yukawa_params[2] = {g, mu};
  double k = 1.0;

  printf("   Born approximation, Yukawa potential (g=%.1f, mu=%.1f, k=%.1f)\n",
         g, mu, k);
  printf("   \\theta   f(\\theta)        d\\sigma/d\\Omega\n");
  printf("   -----   -------------   -------------\n");
  for (double theta = 0.2; theta < M_PI; theta += 0.4) {
    complex_t f =
        born_amplitude(V_yukawa, yukawa_params, k, theta, 50.0, 5000, HBAR_2M);
    double sigma = born_cross_section(f);
    printf("   %5.2f   %13.6e   %13.6e\n", theta, f.re, sigma);
  }

  return 0;
}
