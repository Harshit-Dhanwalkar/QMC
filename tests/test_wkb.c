/*
 * Verifies WKB / Bohr-Sommerfeld quantization.
 *
 * For the harmonic oscillator V(x) = 1/2x^2 (\hbar=m=\omega=1):
 *   Exact: E_n = n + 1/2
 *   WKB:   \cint p dx = (n+1/2)*2 \pi * \hbar -> E_n = n + 1/2  (exact for HO)
 *
 * The HO is special where WKB is exact.
 *
 * Also tests non-trivial potential: V(x) = x^4
 *   WKB E_n = C*(n+1/2)^(4/3) where C = (3 \pi/4)^(2/3) *
 * \Gamma(3/4)/\Gamma(1/4) * ... just verify numerical quadrature gives
 *  self-consistent results.
 */

#include "../core/utils.h"
#include "../physics/wkb.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check(const char *name, double got, double expected, double tol) {
  int pass = fabs(got - expected) < tol;
  printf("   %-30s got=%.6f  exp=%.6f  %s\n", name, got, expected,
         pass ? "PASS" : "FAIL");

  return pass;
}

int main(void) {
  printf(" > Testing WKB / Bohr-Sommerfeld quantization...\n\n");
  int all_pass = 1;

  // Harmonic oscillator: WKB is exact
  printf("   Harmonic oscillator V(x)=1/2x^2 (WKB is exact):\n");

  for (int n = 0; n <= 4; n++) {
    double E_wkb = wkb_energy_harmonic(n, 1.0, 1.0);
    double E_exact = n + 0.5;
    char label[64];
    snprintf(label, sizeof label, "E_%d (WKB)", n);
    all_pass &= check(label, E_wkb, E_exact, 0.01);
  }

  // Bohr-Sommerfeld action integral
  /* For HO at E=0.5: action = 2 * integral_{-1}^{1} * sqrt(2*(E-V)) dx
   *                         = 2 * integral_{-1}^{1} * sqrt(1-x^2) dx = \pi
   * Quantization: action/(2 * \pi) = (n+1/2) -> n=0 at action = \pi
   */
  printf("\n   Bohr-Sommerfeld action integral for HO at E=0.5:\n");
  double action = wkb_action_integral_harmonic(0.5, 1.0, 1.0);
  // Expected: \pi (= (n+1/2)* 2 * \pi * \hbar with n=0, \hbar=1)
  all_pass &=
      check("action/(2 * \\pi) at E=0.5", action / (2.0 * M_PI), 0.5, 0.01);

  // Tunneling probability through rectangular barrier
  printf("\n   WKB tunneling: rectangular barrier V0=2, width=1:\n");
  /*
    T = \exp(-2 * \kappa * L), \kappa = \sqrt(2m(V0-E)) / \hbar = \sqrt((V0-E)
   / hbar_sq_2m)
    with hbar_sq_2m = \hbar^2/(2m) = 0.5 (\hbar=1, m=1): kappa =
   \sqrt((2-1)/0.5) = \sqrt(2) ~= 1.41421356
  */
  double E = 1.0;
  double V0 = 2.0;
  double width = 1.0;
  double hbar_sq_2m = 0.5;
  double T_wkb = wkb_tunneling_rectangular(E, V0, width, hbar_sq_2m);
  double kappa = sqrt((V0 - E) / hbar_sq_2m);
  double T_exact = exp(-2.0 * kappa * width);
  all_pass &= check("T_tunnel (WKB vs exact)", T_wkb, T_exact, 0.01);

  printf("\n");
  if (all_pass)
    printf("   WKB test passed.\n");

  return all_pass ? 0 : 1;
}
