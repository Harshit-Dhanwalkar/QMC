/*
 * Test: Variational Monte Carlo for helium ground state.
 *
 * 1. vmc_trial_wavefunction / vmc_local_energy at fixed configurations must
 *    match reference values -> these are deterministic, no Monte Carlo noise
 *    involved.
 * 2. vmc_run at a reasonable fixed b must satisfy variational theorem (E >=
 *    -2.9037 Hartree, experimental ground state) and fall below plain
 *    product-orbital result (-2.84765625 Hartree, from
 *    helium_ground_state_energy_analytic(2.0)) since Jastrow factor recovers
 *    some correlation energy.
 * 3. Acceptance rates from vmc_run must be in a sane range (order 30-60%),
 *    confirming step sizes are reasonably tuned and sampler isn't degenerate
 *    (always/never accepting).
 * 4. vmc_optimize_b must find a minimum energy at or below the b=0 (pure
 *    product-orbital) energy, with b_opt in physically expected range for this
 *    ansatz.
 */

#include "../physics/helium.h"
#include "../physics/vmc.h"
#include <math.h>
#include <stdio.h>

#ifndef RUNNING_ON_VALGRIND
#define RUNNING_ON_VALGRIND 0
#endif

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

// Deterministic fixture values
static void test_local_energy_fixtures(void) {
  printf("test_local_energy_fixtures:\n");

  {
    vmc_walker_t w = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
    double b = 0.2;
    double E_L = vmc_local_energy(&w, 2.0, 2.0, b);
    double psi = vmc_trial_wavefunction(&w, 2.0, b);
    check_close(E_L, -2.7268846314035815, 1e-9, "E_L case 1");
    check_close(psi, 0.03741385136723659, 1e-9, "Psi case 1");
  }
  {
    vmc_walker_t w = {{0.5, 0.3, -0.2}, {-0.4, 0.6, 0.1}};
    double b = 0.15;
    double E_L = vmc_local_energy(&w, 2.0, 2.0, b);
    double psi = vmc_trial_wavefunction(&w, 2.0, b);
    check_close(E_L, -2.6861570507734376, 1e-9, "E_L case 2");
    check_close(psi, 0.10476678247491533, 1e-9, "Psi case 2");
  }
  {
    vmc_walker_t w = {{2.0, 1.0, 0.5}, {0.3, -0.8, 1.2}};
    double b = 0.35;
    double E_L = vmc_local_energy(&w, 2.0, 2.0, b);
    double psi = vmc_trial_wavefunction(&w, 2.0, b);
    check_close(E_L, -3.3300641407479987, 1e-9, "E_L case 3");
    check_close(psi, 0.0010574877808528424, 1e-9, "Psi case 3");
  }
}

// Degenerate configuration guard (r1 -> 0 must not produce NaN/\infty)
static void test_degenerate_guard(void) {
  printf("test_degenerate_guard:\n");

  vmc_walker_t w = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
  double E_L = vmc_local_energy(&w, 2.0, 2.0, 0.2);

  check_true(isfinite(E_L) && E_L == 0.0,
             "r1->0 configuration returns finite guard value (0.0)");
}

static void test_vmc_run_bounds(void) {
  printf("test_vmc_run_bounds:\n");

  double Z = 2.0, Zeff = 2.0;
  double b = 0.15;
  int burn = RUNNING_ON_VALGRIND ? 500 : 2000;
  int samples = RUNNING_ON_VALGRIND ? 10000 : 200000;
  int block_size = RUNNING_ON_VALGRIND ? 100 : 200;

  vmc_result_t r =
      vmc_run(Z, Zeff, b, burn, samples, block_size, 0.9, 0.9, 20260729ULL);

  printf("  E = %.6f +- %.6f Hartree (n_samples=%d)\n", r.mean, r.error,
         r.n_samples);
  printf("  acceptance: electron1=%.3f electron2=%.3f\n", r.acceptance_rate1,
         r.acceptance_rate2);

  double E_experimental = -2.9037;
  double E_product_orbital = helium_ground_state_energy_analytic(2.0);

  check_true(
      r.mean >= E_experimental - 3.0 * r.error,
      "variational theorem: E >= experimental ground state (within 3\\sigma)");
  check_true(
      r.mean < E_product_orbital,
      "Jastrow-correlated energy improves on plain product-orbital result");
  check_true(r.acceptance_rate1 > 0.15 && r.acceptance_rate1 < 0.85,
             "electron-1 acceptance rate in sane range");
  check_true(r.acceptance_rate2 > 0.15 && r.acceptance_rate2 < 0.85,
             "electron-2 acceptance rate in sane range");
  check_true(r.n_samples == samples, "n_samples matches request");
}

static void test_optimize_b(void) {
  printf("test_optimize_b:\n");

  double Z = 2.0, Zeff = 2.0;
  double b_opt;
  int burn = RUNNING_ON_VALGRIND ? 500 : 1000;
  int samples = RUNNING_ON_VALGRIND ? 20000 : 500000;
  int max_iter = RUNNING_ON_VALGRIND ? 100 : 1000;

  double E_opt = vmc_optimize_b(Z, Zeff, 0.0, 0.6, max_iter, samples, 0.9, 0.9,
                                4242ULL, 1e-3, &b_opt);

  vmc_result_t r_b0 =
      vmc_run(Z, Zeff, 0.0, burn, samples, 200, 0.9, 0.9, 4242ULL);

  printf("  b_opt = %.4f, E(b_opt) = %.6f Hartree\n", b_opt, E_opt);
  printf("  E(b=0) = %.6f Hartree (for comparison)\n", r_b0.mean);

  check_true(b_opt > 0.0 && b_opt < 0.5,
             "optimized b falls in the physically expected range (0, 0.5)");
  check_true(E_opt <= r_b0.mean,
             "optimized energy is at or below the b=0 (no-Jastrow) energy");
  check_true(E_opt >= -2.9037 - 0.02,
             "optimized energy respects the variational bound (with tolerance "
             "for MC noise)");
}

static void test_local_energy_z_neq_zeff_fixture(void) {
  printf("test_local_energy_z_neq_zeff_fixture:\n");

  vmc_walker_t w = {{0.5, 0.3, -0.2}, {-0.4, 0.6, 0.1}};
  double Z = 1.0, Zeff = 2.0; // H- (Z=1) with Zeff=2 trial orbital
  double b = 0.10;
  double E_L = vmc_local_energy(&w, Z, Zeff, b);

  check_close(E_L, 0.29101849950066816, 1e-6, "E_L, Z=1 Zeff=2 (H-)");
}

// Check VMC for a different two-electron ion (Li+, Z=3) lands in a sensible
// range, respecting variational theorem against known-precise reference
// (Reference: -7.2799133 Hartree, Frolov 1997 SVM calculation) with a tolerance
// appropriate for simple Slater-Jastrow
static void test_vmc_run_different_ion(void) {
  printf("test_vmc_run_different_ion:\n");

  double Z = 3.0, Zeff = 2.6; // Li+
  double b = 0.12;
  int burn = RUNNING_ON_VALGRIND ? 500 : 2000;
  int samples = RUNNING_ON_VALGRIND ? 10000 : 20000;
  int block_size = RUNNING_ON_VALGRIND ? 100 : 200;

  vmc_result_t r =
      vmc_run(Z, Zeff, b, burn, samples, block_size, 0.9, 0.9, 99110ULL);

  printf("  Li+ (Z=3, Zeff=%.2f): E = %.6f +- %.6f Hartree\n", Zeff, r.mean,
         r.error);

  double E_exact_liplus = -7.2799133;
  check_true(r.mean >= E_exact_liplus - 3.0 * r.error,
             "Li+ variational theorem: E >= exact ground state (within "
             "3\\sigma)");
  check_true(r.mean < 0.0, "Li+ energy is a sane bound state (negative)");
}

// n_replicas=1 must reproduce vmc_run exactly: vmc_run_parallel's stream 0 is
// rng_seed(master_seed) with zero rng_jump() calls applied, identical to what
// vmc_run does internally, and both route through same vmc_run_with_rng()
// sampling loop.
static void test_vmc_run_parallel_matches_serial_at_one_replica(void) {
  printf("test_vmc_run_parallel_matches_serial_at_one_replica:\n");

  double Z = 2.0, Zeff = 1.6875;
  double b = 0.3;
  int burn = RUNNING_ON_VALGRIND ? 500 : 1000;
  int samples = RUNNING_ON_VALGRIND ? 20000 : 500000;
  int block_size = RUNNING_ON_VALGRIND ? 100 : 200;

  vmc_result_t serial =
      vmc_run(Z, Zeff, b, burn, samples, block_size, 0.9, 0.9, 42ULL);
  vmc_result_t parallel = vmc_run_parallel(1, Z, Zeff, b, burn, samples,
                                           block_size, 0.9, 0.9, 42ULL);

  check_close(parallel.mean, serial.mean, 1e-12,
              "n_replicas=1 mean bit-identical to vmc_run");
  check_close(parallel.acceptance_rate1, serial.acceptance_rate1, 1e-12,
              "n_replicas=1 acceptance_rate1 bit-identical to vmc_run");
}

// Multiple replicas must be genuinely independent samples (not silently same
// chain replayed): total n_samples should be n_replicas * n_samples, and
// combined mean should still satisfy He's variational theorem.
static void test_vmc_run_parallel_helium(void) {
  printf("test_vmc_run_parallel_helium:\n");

  double Z = 2.0, Zeff = 1.6875;
  double b = 0.3;
  int n_replicas = RUNNING_ON_VALGRIND ? 2 : 8;
  int n_samples = RUNNING_ON_VALGRIND ? 20000 : 40000;
  int burn = RUNNING_ON_VALGRIND ? 500 : 1000;
  int block_size = RUNNING_ON_VALGRIND ? 100 : 200;

  vmc_result_t r = vmc_run_parallel(n_replicas, Z, Zeff, b, burn, n_samples,
                                    block_size, 0.9, 0.9, 777ULL);

  printf(
      "  parallel He (n_replicas=%d): E = %.6f +- %.6f Hartree, n_samples=%d\n",
      n_replicas, r.mean, r.error, r.n_samples);

  check_true(r.n_samples == n_replicas * n_samples,
             "n_samples = n_replicas * n_samples");

  double E_exact_helium = -2.9037;
  check_true(r.mean >= E_exact_helium - 3.0 * r.error,
             "parallel He variational theorem: E >= exact ground state (within "
             "3\\sigma)");
  check_true(r.error > 0.0 && r.error < 0.05,
             "parallel He inter-replica error is a sane, nonzero magnitude");
}

int main(void) {
  test_local_energy_fixtures();
  test_local_energy_z_neq_zeff_fixture();
  test_degenerate_guard();
  test_vmc_run_bounds();
  test_vmc_run_different_ion();
  test_optimize_b();
  test_vmc_run_parallel_matches_serial_at_one_replica();
  test_vmc_run_parallel_helium();

  if (failures == 0) {
    printf("\nAll test_vmc checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
