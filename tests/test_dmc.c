/*
 * Test: Diffusion Monte Carlo for helium ground state.
 *
 * 1. dmc_drift_velocity at fixed configurations (deterministic, no MC noise).
 * 2. Full dmc_run must land close to exact helium ground state energy (-2.9037
 *    Hartree) - should be much closer than VMC's variational estimate
 *    (-2.84765625 from plain product orbital, or ~-2.878 from VMC with same
 *    Jastrow).
 * 3. Acceptance rate and mean population must be in a sane range, confirming
 *    sampler isn't degenerate and population control is holding near target.
 */

#include "../core/random.h"
#include "../physics/dmc.h"
#include "../physics/helium.h"
#include "../physics/vmc.h"
#include <math.h>
#include <omp.h>
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

static void test_drift_velocity_fixtures(void) {
  printf("test_drift_velocity_fixtures:\n");

  {
    vmc_walker_t w = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
    double b = 0.2;
    double d0[3], d1[3];

    dmc_drift_params_t p0 = {&w, 0, 2.0, b};
    dmc_drift_params_t p1 = {&w, 1, 2.0, b};

    dmc_drift_velocity(&p0, d0);
    dmc_drift_velocity(&p1, d1);
    check_close(d0[0], -1.7448979591836733, 1e-9, "drift0[0] case 1");
    check_close(d0[1], 0.0, 1e-9, "drift0[1] case 1");
    check_close(d1[0], 1.7448979591836733, 1e-9, "drift1[0] case 1");
  }
  {
    vmc_walker_t w = {{0.5, 0.3, -0.2}, {-0.4, 0.6, 0.1}};
    double b = 0.15;
    double d0[3], d1[3];

    dmc_drift_params_t p0 = {&w, 0, 2.0, b};
    dmc_drift_params_t p1 = {&w, 1, 2.0, b};

    dmc_drift_velocity(&p0, d0);
    dmc_drift_velocity(&p1, d1);
    check_close(d0[0], -1.2797877515194125, 1e-9, "drift0[0] case 2");
    check_close(d0[1], -1.0874706800473128, 1e-9, "drift0[1] case 2");
    check_close(d0[2], 0.5347435312603125, 1e-9, "drift0[2] case 2");
    check_close(d1[0], 0.7564580518012993, 1e-9, "drift1[0] case 2");
    check_close(d1[1], -1.5341846141215307, 1e-9, "drift1[1] case 2");
    check_close(d1[2], -0.16057897463464044, 1e-9, "drift1[2] case 2");
  }
}

static void test_degenerate_guard(void) {
  printf("test_degenerate_guard:\n");

  vmc_walker_t w = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
  double d[3];

  dmc_drift_params_t p0 = {&w, 0, 2.0, 0.2};

  dmc_drift_velocity(&p0, d);
  check_true(isfinite(d[0]) && isfinite(d[1]) && isfinite(d[2]) &&
                 d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.0,
             "r1->0 configuration returns finite zero-vector guard");
}

static void test_dmc_run_accuracy(void) {
  printf("test_dmc_run_accuracy:\n");

  double Z = 2.0, Zeff = 2.0;
  double b = 0.15;
  double tau = 0.01;
  int target_population = RUNNING_ON_VALGRIND ? 20 : 200;
  int max_population = RUNNING_ON_VALGRIND ? 60 : 600;
  int n_equilibration = RUNNING_ON_VALGRIND ? 100 : 500;
  int n_blocks = RUNNING_ON_VALGRIND ? 2 : 20;
  int block_size = RUNNING_ON_VALGRIND ? 20 : 200;

  dmc_result_t r = dmc_run(Z, Zeff, b, target_population, max_population, tau,
                           n_equilibration, n_blocks, block_size,
                           /*seed=*/20260731ULL);

  printf("  mixed:  %.6f +- %.6f Hartree\n", r.energy_mixed, r.error_mixed);
  printf("  growth: %.6f +- %.6f Hartree\n", r.energy_growth, r.error_growth);
  printf("  mean_population=%.2f  acceptance=%.4f  n_blocks=%d\n",
         r.mean_population, r.acceptance_rate, r.n_blocks);

  double E_exact = -2.9037;
  double E_product_orbital = helium_ground_state_energy_analytic(2.0);

  /* NOTE: DMC should land close to exact (much tighter than a VMC-level
   * tolerance would be, but still allowing for residual
   * timestep/population-control bias at finite \tau and finite population). */
  // HACK: Relax tolerance under Valgrind (small population → larger statistical
  // error)
  double tol = RUNNING_ON_VALGRIND ? 0.1 : 0.02;

  check_close(r.energy_mixed, E_exact, tol,
              "mixed estimator close to exact ground state");

  // Skip improvments on VMC checks under Valgrind (too noisy)
  if (!RUNNING_ON_VALGRIND) {
    double E_vmc_jastrow = -2.878; // aprrox value from VMC b-scan
    check_true(r.energy_mixed < E_vmc_jastrow,
               "DMC energy improves on VMC-with-Jastrow estimate");
    check_true(r.energy_mixed < E_product_orbital,
               "DMC energy improves on plain product-orbital estimate");
  }

  check_true(r.energy_mixed > E_exact - 0.05,
             "DMC energy does not undershoot exact value implausibly");

  check_true(r.acceptance_rate > 0.7 && r.acceptance_rate <= 1.0,
             "acceptance rate in reasonable range (expect high, ~0.9+, for "
             "Metropolis-corrected drift-diffusion at this \\tau)");

  check_true(r.n_blocks == n_blocks, "blocks completed");
  check_true(fabs(r.mean_population - target_population) <
                 target_population * 0.1,
             "population control holds near target");
}

/*
 * DMC for a different two-electron ion (Li+, Z=3): must output closer to exact
 * reference than VMC's own Li+ estimate does, and must respect variational
 * theorem.
 */
static void test_dmc_run_different_ion(void) {
  printf("test_dmc_run_different_ion:\n");

  double Z = 3.0, Zeff = 2.6; // Li+
  double b = 0.12;
  double tau = 0.01;
  int target_population = RUNNING_ON_VALGRIND ? 100 : 1000;
  int max_population = RUNNING_ON_VALGRIND ? 1000 : 10000;
  int n_equilibration = RUNNING_ON_VALGRIND ? 100 : 500;
  int n_blocks = RUNNING_ON_VALGRIND ? 2 : 20;
  int block_size = RUNNING_ON_VALGRIND ? 20 : 200;

  // VMC parameters
  int vmc_burn = RUNNING_ON_VALGRIND ? 500 : 2000;
  int vmc_samples = RUNNING_ON_VALGRIND ? 10000 : 100000;
  int vmc_block = 100;

  double E_exact_liplus = -7.2799133;

  vmc_result_t vmc_r =
      vmc_run(Z, Zeff, b, vmc_burn, vmc_samples, vmc_block, 0.9, 0.9, 55221ULL);
  dmc_result_t dmc_r = dmc_run(Z, Zeff, b, target_population, max_population,
                               tau, n_equilibration, n_blocks, block_size,
                               /*seed=*/55222ULL);

  printf("  Li+ VMC: %.6f +- %.6f Hartree\n", vmc_r.mean, vmc_r.error);
  printf("  Li+ DMC (mixed): %.6f +- %.6f Hartree\n", dmc_r.energy_mixed,
         dmc_r.error_mixed);

  check_true(dmc_r.energy_mixed >= E_exact_liplus - 3.0 * dmc_r.error_mixed,
             "Li+ DMC respects the variational theorem");

  // HACK: Under Valgrind, skip the improvement check (too noisy)
  if (!RUNNING_ON_VALGRIND) {
    check_true(dmc_r.energy_mixed < vmc_r.mean,
               "Li+ DMC improves on Li+ VMC estimate (projects toward exact)");
  }
}

/*
 * Regression test for max_population set close enough to target_population that
 * population control should trigger "on nearly every generation"): with
 * max_population only 3x target_population and a small tau=0.01, the E_T
 * feedback controller alone keeps the population close enough to target that
 * comb resampling is essentially never needed.
 */
static void test_dmc_frequent_resampling(void) {
  printf("test_dmc_frequent_resampling:\n");

  double Z = 2.0, Zeff = 2.0;
  double b = 0.15;
  double tau = 0.01;
  int target_population = RUNNING_ON_VALGRIND ? 20 : 200;
  int max_population = RUNNING_ON_VALGRIND ? 60 : 600;
  int n_equilibration = RUNNING_ON_VALGRIND ? 100 : 500;
  int n_blocks = RUNNING_ON_VALGRIND ? 2 : 20;
  int block_size = RUNNING_ON_VALGRIND ? 20 : 200;

  dmc_result_t r = dmc_run(Z, Zeff, b, target_population, max_population, tau,
                           n_equilibration, n_blocks, block_size,
                           /*seed=*/20260803ULL);

  printf("  mixed:  %.6f +- %.6f Hartree  (max_population=%d, target=%d)\n",
         r.energy_mixed, r.error_mixed, max_population, target_population);
  printf("  mean_population=%.2f  acceptance=%.4f  n_resamples=%ld\n",
         r.mean_population, r.acceptance_rate, r.n_resamples);

  double E_exact = -2.9037;
  double tol = RUNNING_ON_VALGRIND ? 0.1 : 0.03;
  check_close(r.energy_mixed, E_exact, tol,
              "DMC still lands close to exact under frequent resampling");
  check_true(fabs(r.mean_population - target_population) <
                 target_population * 0.1,
             "population control holds near target under frequent triggering");
}

/*
 * Stress test for comb-resampling code path itself (not just E_T feedback
 * controller): forces post-branching population overflow past max_population by
 * combining a tight target/max_population margin with a large \tau (unphysical
 * for accurate DMC energies : this is purely a population-control mechanism
 * test, not a physics accuracy test).
 */
static void test_dmc_resampling_engages_under_tight_margin(void) {
  printf("test_dmc_resampling_engages_under_tight_margin:\n");

  double Z = 2.0, Zeff = 1.6875;
  double b = 0.35;
  double tau = 0.1; /* deliberately large: amplifies branching-weight
                     * variance and weakens the kappa/tau feedback gain, so
                     * population pressure against the cap is easy to
                     * reproduce deterministically for a fixed seed */
  int target_population = RUNNING_ON_VALGRIND ? 30 : 150;
  int max_population =
      RUNNING_ON_VALGRIND
          ? 33
          : 165; /* only 1.1x target, far below documented-safe ~3x */
  int n_equilibration = RUNNING_ON_VALGRIND ? 20 : 100;
  int n_blocks = RUNNING_ON_VALGRIND ? 3 : 10;
  int block_size = RUNNING_ON_VALGRIND ? 20 : 100;

  dmc_result_t r = dmc_run(Z, Zeff, b, target_population, max_population, tau,
                           n_equilibration, n_blocks, block_size,
                           /*seed=*/7ULL);

  printf("  n_resamples=%ld  mean_population=%.2f  energy_mixed=%.4f\n",
         r.n_resamples, r.mean_population, r.energy_mixed);

  check_true(
      r.n_resamples > 0,
      "comb resampling actually engages under a tight max_population margin");
  check_true(r.mean_population > target_population * 0.1 &&
                 r.mean_population < max_population * 1.5,
             "population stays bounded (neither collapsed to ~0 nor "
             "unboundedly growing) with resampling active");
  check_true(r.energy_mixed == r.energy_mixed,
             "energy_mixed is not NaN (sanity check on stressed run)");
}

/*
 * n_replicas=1 must reproduce dmc_run exactly, same reason being as analogous
 * VMC regression test: stream 0 has zero rng_jump() calls applied, and both
 * route through same dmc_run_with_rng() population loop. (same parameters, same
 * seed)
 */
static void test_dmc_run_parallel_matches_serial_at_one_replica(void) {
  printf("test_dmc_run_parallel_matches_serial_at_one_replica:\n");

  double Z = 2.0, Zeff = 2.0;
  double b = 0.15;
  double tau = 0.01;
  int target_population = RUNNING_ON_VALGRIND ? 20 : 200;
  int max_population = RUNNING_ON_VALGRIND ? 60 : 600;
  int n_equilibration = RUNNING_ON_VALGRIND ? 100 : 200;
  int n_blocks = RUNNING_ON_VALGRIND ? 5 : 20;
  int block_size = RUNNING_ON_VALGRIND ? 20 : 200;

  dmc_result_t serial =
      dmc_run(Z, Zeff, b, target_population, max_population, tau,
              n_equilibration, n_blocks, block_size, /*seed=*/55ULL);
  dmc_result_t parallel =
      dmc_run_parallel(1, Z, Zeff, b, target_population, max_population, tau,
                       n_equilibration, n_blocks, block_size, /*seed=*/55ULL);

  check_close(parallel.energy_mixed, serial.energy_mixed, 1e-12,
              "n_replicas=1 energy_mixed bit-identical to dmc_run");
  check_close(parallel.mean_population, serial.mean_population, 1e-12,
              "n_replicas=1 mean_population bit-identical to dmc_run");
}

/*
 * NOTE: The intra-run walker population loop is OpenMP-parallelized across
 * walkers within a single DMC run (not just across independent replicas, as
 * dmc_run_parallel already did). Each population slot `i` always draws from its
 * own fixed walker_streams[i] (derived once via rng_jump chaining before
 * parallel region opens, independent of which thread ends up processing that
 * lot in any given generation), so results should be 'bit-identical' regardless
 * of how many threads OpenMP actually uses : a correctness property for
 * thread-safety of this refactor (not just "close within statistical noise").
 */
static void test_dmc_run_deterministic_across_thread_counts(void) {
  printf("test_dmc_run_deterministic_across_thread_counts:\n");

  double Z = 2.0, Zeff = 2.0;
  double b = 0.15;
  double tau = 0.01;
  int target_population = RUNNING_ON_VALGRIND ? 20 : 150;
  int max_population = RUNNING_ON_VALGRIND ? 60 : 450;
  int n_equilibration = RUNNING_ON_VALGRIND ? 50 : 100;
  int n_blocks = RUNNING_ON_VALGRIND ? 5 : 10;
  int block_size = RUNNING_ON_VALGRIND ? 20 : 100;

  int max_threads = omp_get_max_threads();

  omp_set_num_threads(1);
  dmc_result_t one_thread =
      dmc_run(Z, Zeff, b, target_population, max_population, tau,
              n_equilibration, n_blocks, block_size, /*seed=*/777ULL);

  omp_set_num_threads(max_threads > 1 ? max_threads : 2);
  dmc_result_t many_threads =
      dmc_run(Z, Zeff, b, target_population, max_population, tau,
              n_equilibration, n_blocks, block_size, /*seed=*/777ULL);

  omp_set_num_threads(max_threads); /* restore for any later tests */

  printf("  max_threads available=%d\n", max_threads);
  check_close(many_threads.energy_mixed, one_thread.energy_mixed, 1e-12,
              "energy_mixed bit-identical across thread counts");
  check_close(many_threads.mean_population, one_thread.mean_population, 1e-12,
              "mean_population bit-identical across thread counts");
  check_close(many_threads.acceptance_rate, one_thread.acceptance_rate, 1e-12,
              "acceptance_rate bit-identical across thread counts");
}

/*
 * Parallel DMC with multiple replicas
 * Multiple independent DMC populations combined must still land close to the
 * exact helium ground state, with a sane nonzero inter-replica error bar.
 */
static void test_dmc_run_parallel_helium(void) {
  printf("test_dmc_run_parallel_helium:\n");

  int n_replicas = RUNNING_ON_VALGRIND ? 2 : 6;
  double Z = 2.0, Zeff = 2.0;
  double b = 0.15;
  double tau = 0.01;
  int target_population = RUNNING_ON_VALGRIND ? 20 : 200;
  int max_population = RUNNING_ON_VALGRIND ? 60 : 600;
  int n_blocks = RUNNING_ON_VALGRIND ? 2 : 20;
  int block_size = RUNNING_ON_VALGRIND ? 20 : 200;
  int n_equilibration = RUNNING_ON_VALGRIND ? 100 : 500;

  dmc_result_t r = dmc_run_parallel(n_replicas, Z, Zeff, b, target_population,
                                    max_population, tau, n_equilibration,
                                    n_blocks, block_size,
                                    /*master_seed=*/31415ULL);

  printf("  parallel He (n_replicas=%d): mixed=%.6f +- %.6f Hartree\n",
         n_replicas, r.energy_mixed, r.error_mixed);

  double E_exact = -2.9037;

  double tol = RUNNING_ON_VALGRIND ? 0.1 : 0.03;
  check_close(r.energy_mixed, E_exact, tol,
              "parallel DMC lands close to exact helium ground state");

  check_true(r.error_mixed > 0.0 &&
                 r.error_mixed < (RUNNING_ON_VALGRIND ? 0.1 : 0.05),
             "parallel DMC inter-replica error is a sane, nonzero magnitude");
  check_true(r.n_blocks == n_blocks * n_replicas,
             "n_blocks reports n_blocks_per_replica * n_replicas");
}

int main(void) {
  test_drift_velocity_fixtures();
  test_degenerate_guard();
  test_dmc_run_accuracy();
  test_dmc_frequent_resampling();
  test_dmc_resampling_engages_under_tight_margin();
  test_dmc_run_different_ion();
  test_dmc_run_parallel_matches_serial_at_one_replica();
  test_dmc_run_deterministic_across_thread_counts();
  test_dmc_run_parallel_helium();

  if (failures == 0) {
    printf("\nAll test_dmc checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
