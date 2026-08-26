/*
 * Test: physics/wavefn.c's momentum-space expectation values
 * (wavefunction_expect_p, wavefunction_expect_p2, wavefunction_delta_p),
 * and physics/uncertainty.c's compute_uncertainties/compute_energy_expectation.
 *
 *   1. core/utils.c's position_to_momentum called fft_shift() after the
 *      FFT, but its only two consumers (wavefunction_expect_p/p2 below)
 *      both index the result using the standard unshifted fft bin order.
 *   2. core/utils.c's position_to_momentum used fft_normalized's unitary-DFT (1
 *      / \sqrt(N)) scaling convention directly, with no further rescaling
 *      position-to-momentum-space normalization for a continuous wavefunction
 *      needs an explicit dx-dependent prefactor (\phi(k) = (dx / \sqrt(2 *
 *      \pi)) * FFT_raw(\psi_x)), "signal processing" convention.
 */

#include "../physics/uncertainty.h"
#include "../physics/wavefn.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  if (fabs(got - expected) > tol) {
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

/* Build the exact QHO ground state (natural units: hbar=m=omega=1),
 * \psi(x) = \pi^(-1/4) * \exp(-x^2/2), on a grid of n points over [-L/2,L/2].
 */
static wavefunction_t *make_qho_ground_state(int n, double L) {
  wavefunction_t *wf = wavefunction_alloc(n);
  wf->dx = L / (n - 1);

  double norm_const = pow(1.0 / M_PI, 0.25);
  for (int i = 0; i < n; i++) {
    double xi = -L / 2.0 + i * wf->dx;

    wf->x[i] = xi;
    wf->psi->data[i].re = norm_const * exp(-xi * xi / 2.0);
    wf->psi->data[i].im = 0.0;
  }

  return wf;
}

static void test_qho_momentum_expectations_power_of_two_grid(void) {
  printf("Test: QHO ground state momentum expectations on a power-of-two grid "
         "(radix-2 FFT path)\n");

  wavefunction_t *wf = make_qho_ground_state(2048, 16.0);

  double norm = 0.0;
  for (int i = 0; i < wf->n; i++) {
    norm += wf->psi->data[i].re * wf->psi->data[i].re * wf->dx;
  }

  check_close(norm, 1.0, 1e-6, "position-space normalization holds");

  check_close(wavefunction_expect_p(wf), 0.0, 1e-6,
              "<p>=0 for the real, even ground state");
  check_close(wavefunction_expect_p2(wf), 0.5, 1e-3,
              "<p^2>=0.5 exactly (natural units, \\hbar=m=\\omega=1)");
  check_close(wavefunction_delta_p(wf), sqrt(0.5), 1e-3,
              "\\Delta_p = \\sqrt(<p^2>-<p>^2) = \\sqrt(0.5)");

  wavefunction_free(wf);
}

static void test_qho_momentum_expectations_non_power_of_two_grid(void) {
  printf("Test: QHO ground state momentum expectations on a NON-power-of-two "
         "grid\n");

  wavefunction_t *wf = make_qho_ground_state(2001, 16.0);

  check_close(wavefunction_expect_p2(wf), 0.5, 1e-2,
              "<p^2>=0.5 on a non-power-of-two grid too");

  wavefunction_free(wf);
}

static void test_uncertainty_principle_saturated(void) {
  printf("Test: compute_uncertainties on the QHO ground state saturates "
         "Heisenberg bound exactly, \\Delta_x * \\Delta_p = \\hbar/2 = 0.5\n");

  wavefunction_t *wf = make_qho_ground_state(2048, 16.0);
  uncertainty_t u = compute_uncertainties(wf);

  check_close(u.mean_x, 0.0, 1e-6, "<x>=0");
  check_close(u.mean_p, 0.0, 1e-6, "<p>=0");
  check_close(u.delta_x, sqrt(0.5), 1e-6, "\\Delta_x = \\sqrt(0.5)");
  check_close(u.delta_p, sqrt(0.5), 1e-3, "\\Delta_p = \\sqrt(0.5)");
  check_close(
      u.delta_x * u.delta_p, 0.5, 1e-3,
      "\\Delta_x * \\Delta_p = \\hbar/2 exactly (minimum-uncertainty state)");

  wavefunction_free(wf);
}

static double harmonic_V(double x, void *params) {
  (void)params;
  return 0.5 * x * x;
}

static void test_energy_expectation_matches_exact_ground_energy(void) {
  printf("Test: compute_energy_expectation on the QHO ground state matches "
         "exact ground-state energy \\hbar * \\omega/2 = 0.5\n");

  wavefunction_t *wf = make_qho_ground_state(2048, 16.0);

  double H = compute_energy_expectation(wf, harmonic_V, NULL, 1.0);
  check_close(H, 0.5, 1e-3,
              "<H> = <T>+<V> = 0.25+0.25 = 0.5, matching the exact "
              "ground-state energy\n");

  wavefunction_free(wf);
}

int main(void) {
  printf("=== physics/wavefn.c momentum-expectation-value tests ===\n\n");

  test_qho_momentum_expectations_power_of_two_grid();
  test_qho_momentum_expectations_non_power_of_two_grid();
  test_uncertainty_principle_saturated();
  test_energy_expectation_matches_exact_ground_energy();

  if (failures == 0) {
    printf("\nAll test_wavefn checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_wavefn check(s) FAILED.\n", failures);
    return 1;
  }
}
