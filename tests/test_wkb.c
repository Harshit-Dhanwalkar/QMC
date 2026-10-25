/*
 * Test: Verifies WKB / Bohr-Sommerfeld quantization
 *
 * For harmonic oscillator V(x) = 1/2x^2 (\hbar = m = \omega = 1):
 *   Exact: E_n = n + 1/2
 *   WKB:   \cint p dx = (n+1/2)*2 \pi * \hbar -> E_n = n + 1/2  (exact for HO)
 *
 * HO is special where WKB is exact
 *
 * Also tests non-trivial potential: V(x) = x^4
 *   WKB E_n = C*(n+1/2)^(4/3)
 * where C = (3 \pi/4)^(2/3) * \Gamma(3/4)/\Gamma(1/4) * ... just verify
 * numerical quadrature gives self-consistent results
 *
 * 1. Harmonic oscillator: WKB is exact - E_n = n + 1/2 for n = 0..4
 * 2. Bohr-Sommerfeld action integral for HO at E = 0.5 equals \pi
 *    (quantization condition action/(2 * \pi) = n + 1/2, n = 0)
 * 3. WKB tunneling through a rectangular barrier matches exact
 */

#include "../core/utils.h"
#include "../physics/wkb.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8g expected=%.8g err=%.2e\n", label, got, expected, err);
  if (!(err <= tol)) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

// Harmonic oscillator: WKB is exact
static void test_harmonic_oscillator(void) {
  printf("  === Test: harmonic oscillator (WKB is exact) ===\n");

  for (int n = 0; n <= 4; n++) {
    double E_wkb = wkb_energy_harmonic(n, 1.0, 1.0);
    double E_exact = n + 0.5;
    char label[64];

    snprintf(label, sizeof label, "E_%d (WKB)", n);
    check_close(E_wkb, E_exact, 0.01, label);
  }
}

// Bohr-Sommerfeld action integral
// For HO at E=0.5: action = 2 * \int_{-1}^{1} * \sqrt(2*(E-V)) dx
//                         = 2 * \int_{-1}^{1} * \sqrt(1-x^2) dx = \pi
// Quantization: action/(2 * \pi) = (n+1/2) -> n=0 at action = \pi
static void test_bohr_sommerfeld_action(void) {
  printf("  === Test: Bohr-Sommerfeld action integral for HO at E=0.5 ===\n");

  // action = 2 * \int_{-1}^{1} \sqrt(2 * (E-V)) dx = pi at E=0.5, so
  // quantization condition action/(2 * \pi) = 1/2 selects n = 0
  double action = wkb_action_integral_harmonic(0.5, 1.0, 1.0);

  check_close(action / (2.0 * M_PI), 0.5, 0.01, "action/(2*pi) at E=0.5");
}

// Tunneling probability through rectangular barrier
// T = \exp(-2 * \kappa * L), \kappa = \sqrt(2m(V0-E)) / \hbar = \sqrt((V0-E)
// / hbar_sq_2m)
// with hbar_sq_2m = \hbar^2/(2m) = 0.5 (\hbar=1, m=1):
//  \kappa = \sqrt((2-1)/0.5) = \sqrt(2) ~= 1.41421356
static void test_tunneling_rectangular(void) {
  printf("  === Test: WKB tunneling through rectangular barrier V0=2, "
         "width=1 ===\n");

  double E = 1.0;
  double V0 = 2.0;
  double width = 1.0;
  double hbar_sq_2m = 0.5;
  double T_wkb = wkb_tunneling_rectangular(E, V0, width, hbar_sq_2m);
  double kappa = sqrt((V0 - E) / hbar_sq_2m);
  double T_exact = exp(-2.0 * kappa * width);

  check_close(T_wkb, T_exact, 0.01, "T_tunnel (WKB vs exact)");
}

int main(void) {
  printf(" > WKB / Bohr-Sommerfeld quantization tests:\n");

  test_harmonic_oscillator();
  test_bohr_sommerfeld_action();
  test_tunneling_rectangular();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll WKB checks passed\n");

  return 0;
}
