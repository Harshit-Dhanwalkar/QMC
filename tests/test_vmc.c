/*
Test: Variational Monte Carlo for helium ground state.

1. vmc_trial_wavefunction / vmc_local_energy at fixed configurations must match
  independently-computed (Python/sympy-derived, then finite-difference
  cross-checked) reference values -> these are deterministic, no Monte Carlo
  noise involved.
2. vmc_run at a reasonable fixed b must satisfy variational theorem (E >=
  -2.9037 Hartree, experimental ground state) and fall below plain
  product-orbital result (-2.84765625 Hartree, from
  helium_ground_state_energy_analytic(2.0)) since Jastrow factor recovers some
  correlation energy.
3. Acceptance rates from vmc_run must be in a sane range (order 30-60%),
  confirming step sizes are reasonably tuned and sampler isn't degenerate
  (always/never accepting).
4. vmc_optimize_b must find a minimum energy at or below the b=0 (pure
  product-orbital) energy, with b_opt in physically expected range for this
  ansatz.
*/

#include "../core/random.h"
#include "../physics/helium.h"
#include "../physics/vmc.h"
#include <math.h>
#include <stdio.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.10f expected=%.10f err=%.2e\n", label, got, expected,
         err);
  if (err > tol) {
    printf("  FAIL: %s\n", label);

    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

// Deterministic fixture values, computed independently in Python from same
// closed-form local-energy expression (cross-checked there against
// finite-difference Laplacian to ~1e-7).
static void test_local_energy_fixtures(void) {
  printf("test_local_energy_fixtures:\n");

  {
    vmc_walker_t w = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
    double b = 0.2;
    double E_L = vmc_local_energy(&w, 2.0, b);
    double psi = vmc_trial_wavefunction(&w, 2.0, b);
    check_close(E_L, -2.7268846314035815, 1e-9, "E_L case 1");
    check_close(psi, 0.03741385136723659, 1e-9, "Psi case 1");
  }
  {
    vmc_walker_t w = {{0.5, 0.3, -0.2}, {-0.4, 0.6, 0.1}};
    double b = 0.15;
    double E_L = vmc_local_energy(&w, 2.0, b);
    double psi = vmc_trial_wavefunction(&w, 2.0, b);
    check_close(E_L, -2.6861570507734376, 1e-9, "E_L case 2");
    check_close(psi, 0.10476678247491533, 1e-9, "Psi case 2");
  }
  {
    vmc_walker_t w = {{2.0, 1.0, 0.5}, {0.3, -0.8, 1.2}};
    double b = 0.35;
    double E_L = vmc_local_energy(&w, 2.0, b);
    double psi = vmc_trial_wavefunction(&w, 2.0, b);
    check_close(E_L, -3.3300641407479987, 1e-9, "E_L case 3");
    check_close(psi, 0.0010574877808528424, 1e-9, "Psi case 3");
  }
}

// Degenerate configuration guard (r1 -> 0 must not produce NaN/\infty)
static void test_degenerate_guard(void) {
  printf("test_degenerate_guard:\n");

  vmc_walker_t w = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
  double E_L = vmc_local_energy(&w, 2.0, 0.2);
  check_true(isfinite(E_L) && E_L == 0.0,
             "r1->0 configuration returns finite guard value (0.0)");
}

static void test_vmc_run_bounds(void) {
  printf("test_vmc_run_bounds:\n");

  double Zeff = 2.0;
  double b = 0.15;
  vmc_result_t r = vmc_run(Zeff, b, 2000, 200000, 200, 0.9, 0.9, 20260729ULL);

  printf("  E = %.6f +- %.6f Hartree (n_samples=%d)\n", r.mean, r.error,
         r.n_samples);
  printf("  acceptance: electron1=%.3f electron2=%.3f\n", r.acceptance_rate1,
         r.acceptance_rate2);

  double E_experimental = -2.9037;
  double E_product_orbital = helium_ground_state_energy_analytic(2.0);

  check_true(r.mean >= E_experimental - 3.0 * r.error,
             "variational theorem: E >= experimental ground state (within "
             "3 sigma)");
  check_true(r.mean < E_product_orbital,
             "Jastrow-correlated energy improves on plain product-orbital "
             "result");
  check_true(r.acceptance_rate1 > 0.15 && r.acceptance_rate1 < 0.85,
             "electron-1 acceptance rate in sane range");
  check_true(r.acceptance_rate2 > 0.15 && r.acceptance_rate2 < 0.85,
             "electron-2 acceptance rate in sane range");
  check_true(r.n_samples == 200000, "n_samples matches request (block_size "
                                    "divides evenly here)");
}

static void test_optimize_b(void) {
  printf("test_optimize_b:\n");

  double Zeff = 2.0;
  double b_opt;
  double E_opt = vmc_optimize_b(Zeff, 0.0, 0.6, 1000, 50000, 0.9, 0.9, 4242ULL,
                                1e-3, &b_opt);

  vmc_result_t r_b0 = vmc_run(Zeff, 0.0, 1000, 50000, 200, 0.9, 0.9, 4242ULL);

  printf("  b_opt = %.4f, E(b_opt) = %.6f Hartree\n", b_opt, E_opt);
  printf("  E(b=0) = %.6f Hartree (for comparison)\n", r_b0.mean);

  check_true(b_opt > 0.0 && b_opt < 0.5,
             "optimized b falls in the physically expected range (0, 0.5)");
  check_true(E_opt <= r_b0.mean,
             "optimized energy is at or below the b=0 (no-Jastrow) energy");
  check_true(E_opt >= -2.9037 - 0.02,
             "optimized energy respects the variational bound (with "
             "tolerance for MC noise)");
}

int main(void) {
  test_local_energy_fixtures();
  test_degenerate_guard();
  test_vmc_run_bounds();
  test_optimize_b();

  if (failures == 0) {
    printf("\nAll test_vmc checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
