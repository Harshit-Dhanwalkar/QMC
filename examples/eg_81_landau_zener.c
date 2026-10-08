/*
 * Landau-Zener transitions at an avoided crossing
 *
 * NOTE: A two-level system starts in lower adiabatic state and detuning is
 * swept linearly through resonance, H = (1/2)(\Omega \sigma_x + rate t
 * \sigma_z). Probability of ending in upper level (a diabatic jump) is P =
 * \exp(-\pi \Omega^2 / (2 rate)). Prints it for a range of sweep rates next to
 * RK4 integration, then two-pass (Stueckelberg) interference as sweep amplitude
 * is changed
 */

#include "../physics/landau_zener.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  const double omega = 1.0;
  const double amp = 30.0;

  printf(" > Landau-Zener sweep, Omega = %.1f, sweep +/- %.0f\n\n", omega, amp);
  printf("  %8s  %12s  %12s\n", "rate", "integrated", "exp(-pi O^2/2v)");
  for (double rate = 0.1; rate <= 10.01; rate *= 2.0) {
    double v[3];

    lz_lower_state(v, omega, -amp);
    if (lz_sweep(v, omega, rate, amp, 0.01, 1) != 0) {
      return 1;
    }

    printf("  %8.3f  %12.6f  %12.6f\n", rate,
           lz_upper_population(v, omega, amp), lz_probability(omega, rate));
  }

  double p = lz_probability(omega, 1.0);

  printf("\n  Two passes at rate 1 (P_LZ = %.4f), envelope 4P(1-P) = %.4f\n", p,
         4.0 * p * (1.0 - p));
  for (int k = 0; k < 12; k++) {
    double a = 20.0 + 0.1 * k;
    double v[3];

    lz_lower_state(v, omega, -a);
    lz_sweep(v, omega, 1.0, a, 0.01, 2);

    printf("  amp %5.2f  P_upper = %.4f\n", a,
           lz_upper_population(v, omega, -a));
  }

  return 0;
}
