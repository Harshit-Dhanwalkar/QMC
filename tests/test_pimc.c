/*
 * Test: Path Integral Monte Carlo for helium ground state (Kelbg-regularized
 * Coulomb, bisection sampling).
 *
 * 1. kelbg_potential / kelbg_energy_correction at fixed (r, \lambda, q) must
 *    match independently-derived closed-form values - deterministic, no MC
 *    noise.
 * 2. Kelbg potential must be finite at r=0 and must reduce to the bare Coulomb
 *    potential q/r as \lambda -> 0 (\tau -> 0 high-temperature limit).
 * 3. Zero-charge sanity check: q=0 must give exactly V=0 for both functions, at
 *    any r and \lambda (confirms charge scaling isn't accidentally offset).
 * 4. Invalid-input / edge-case handling for allocation and bisection moves.
 * 5. Full pimc_run output close to exact helium ground state energy (-2.9037
 *    Hartree) at validated (P, \tau, level) parameters, with reasonable
 *    acceptance rate.
 *
 * FIX: Thermodynamic energy estimator's variance grows with P (PIMC issue:
 * kinetic term is difference of two large, partially-cancelling quantities), so
 * pushing P much beyond ~512 at fixed statistics needs proportionally more
 * samples for tight comparison, or a virial estimator (which is not
 * implemented) to avoid growing variance altogether.
 * HACK: This test uses P=512 specifically because it's validated and found
 * during development, not largest P that runs.
 */

#include "../core/random.h"
#include "../physics/pimc.h"
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

static void test_kelbg_fixtures(void) {
  printf("test_kelbg_fixtures:\n");

  {
    double r = 0.618;
    double tau = 0.139;
    double lambda = sqrt(tau / 0.5);
    double q = 1.0;
    double v = kelbg_potential(r, lambda, q);
    double vc = kelbg_energy_correction(r, lambda, q);
    check_close(v, 1.5359326977018246, 1e-9, "kelbg_potential fixture 1");
    check_close(vc, 1.3722250340689774, 1e-9,
                "kelbg_energy_correction fixture 1");
  }
  {
    double r = 2.885;
    double tau = 0.365;
    double lambda = sqrt(tau / 0.5);
    double q = -3.0;
    double v = kelbg_potential(r, lambda, q);
    double vc = kelbg_energy_correction(r, lambda, q);
    check_close(v, -1.0398608977874575, 1e-9, "kelbg_potential fixture 2");
    check_close(vc, -1.0398553136873234, 1e-9,
                "kelbg_energy_correction fixture 2");
  }
}

static void test_kelbg_limits(void) {
  printf("test_kelbg_limits:\n");

  // Finite at r=0: analytic limit is (q * \sqrt(\pi) / \lambda)
  double lambda = 0.2, q = -2.0;
  double v0 = kelbg_potential(1e-6, lambda, q);

  check_true(isfinite(v0), "V_Kelbg finite as r->0");
  check_close(
      v0, q * sqrt(M_PI) / lambda, 1e-3,
      "V_Kelbg(r->0) matches analytic limit q * \\sqrt(\\pi) / \\lambda");

  // Reduces to bare Coulomb q/r as \lambda -> 0 (\tau -> 0)
  double r = 0.5;
  double v_small_lambda = kelbg_potential(r, 1e-3, q);

  check_close(v_small_lambda, q / r, 1e-2,
              "V_Kelbg -> bare Coulomb as \\lambda -> 0");

  // Zero charge sanity: q=0 must give exactly 0.
  check_close(kelbg_potential(0.7, 0.3, 0.0), 0.0, 1e-15,
              "kelbg_potential(q=0) == 0 exactly");
  check_close(kelbg_energy_correction(0.7, 0.3, 0.0), 0.0, 1e-15,
              "kelbg_energy_correction(q=0) == 0 exactly");
}

static void test_walker_alloc_and_init(void) {
  printf("test_walker_alloc_and_init:\n");

  check_true(pimc_walker_alloc(0) == NULL, "P=0 allocation rejected");
  check_true(pimc_walker_alloc(-5) == NULL, "negative P allocation rejected");

  pimc_walker_t *w = pimc_walker_alloc(16);
  check_true(w != NULL, "P=16 allocation succeeds");

  rng_state_t rng;
  rng_seed(&rng, 111ULL);
  pimc_walker_init(w, &rng, 2.0);

  int all_finite = 1;
  for (int i = 0; i < w->P; i++) {
    for (int k = 0; k < 3; k++) {
      if (!isfinite(w->r1[i][k]) || !isfinite(w->r2[i][k])) {
        all_finite = 0;
      }
    }
  }

  check_true(all_finite, "initialized walker has all-finite bead positions");

  pimc_walker_free(w);
}

static void test_bisection_invalid_input(void) {
  printf("test_bisection_invalid_input:\n");

  pimc_walker_t *w = pimc_walker_alloc(16);
  rng_state_t rng;
  rng_seed(&rng, 222ULL);
  pimc_walker_init(w, &rng, 2.0);

  check_true(pimc_bisection_move(NULL, 0, 2.0, 0.05, 2, &rng) == 0,
             "NULL walker rejected");
  check_true(pimc_bisection_move(w, 2, 2.0, 0.05, 2, &rng) == 0,
             "invalid which rejected");
  check_true(pimc_bisection_move(w, 0, 2.0, 0.05, 0, &rng) == 0,
             "level=0 rejected");
  check_true(pimc_bisection_move(w, 0, 2.0, 0.05, 5, &rng) == 0,
             "level > log2(P) rejected (P=16, log2=4)");

  pimc_walker_free(w);
}

// Full run at validated (P, \tau, level) parameters
static void test_pimc_run_helium_accuracy(void) {
  printf("test_pimc_run_helium_accuracy:\n");

  double Z = 2.0;
  double beta = 8.0;
  int P = RUNNING_ON_VALGRIND ? 64 : 512;
  int n_blocks = RUNNING_ON_VALGRIND ? 10 : 40;
  int samples_per_block = RUNNING_ON_VALGRIND ? 100 : 300;
  double tau = beta / P;
  int level = 4;

  pimc_result_t r =
      pimc_run(Z, P, tau, level, 400, n_blocks, samples_per_block, 909090ULL);

  printf("  P=%d \\tau=%.5f \\beta=%.2f: E=%.6f +- %.6f  E_virial=%.6f +- %.6f "
         "acceptance=%.4f  n_blocks=%d\n",
         P, tau, P * tau, r.energy, r.error, r.energy_virial, r.error_virial,
         r.acceptance_rate, r.n_blocks);

  double E_exact = -2.9037;

  check_true(r.n_blocks == n_blocks, "requested number of blocks completed");
  check_true(r.acceptance_rate > 0.3,
             "bisection acceptance stays well above near-zero pathology");

  // HACK: Under Valgrind, reduced statistics give large error; relax tolerance
  // heavily.
  double tol = RUNNING_ON_VALGRIND ? 2.0 : 0.4;
  check_close(r.energy, E_exact, tol,
              "PIMC energy lands close to exact helium ground state");
  check_close(r.energy_virial, E_exact, tol,
              "PIMC virial energy lands close to exact helium ground state");
}

/* NOTE: The virial estimator (pimc_virial_estimator) is a different formula
 * from the thermodynamic one (pimc_energy_estimator), derived independently via
 * a coordinate-rescaling argument. Measured on the same sampled configurations,
 * both must agree within their combined statistical error
 */
static void test_virial_matches_thermodynamic(void) {
  printf("test_virial_matches_thermodynamic:\n");

  double Z = 2.0;
  double beta = 8.0;
  int P = RUNNING_ON_VALGRIND ? 32 : 256;
  int n_blocks = RUNNING_ON_VALGRIND ? 4 : 40;
  int samples_per_block = RUNNING_ON_VALGRIND ? 50 : 300;
  double tau = beta / P;
  int level = 4;

  pimc_result_t r =
      pimc_run(Z, P, tau, level, 400, n_blocks, samples_per_block, 20260805ULL);

  double diff = fabs(r.energy - r.energy_virial);
  double combined_err =
      sqrt(r.error * r.error + r.error_virial * r.error_virial);

  printf("  thermodynamic: %.6f +- %.6f\n", r.energy, r.error);
  printf("  virial:        %.6f +- %.6f\n", r.energy_virial, r.error_virial);
  printf("  |diff|=%.6f  combined_err (1\\sigma)=%.6f\n", diff, combined_err);

  check_true(diff < 4.0 * combined_err,
             "thermodynamic and virial estimators agree within their combined "
             "statistical error");
}

/* NOTE: The motivation for virial estimator is lower variance than
 * thermodynamic one at same sample size.
 */
static void test_virial_lower_variance(void) {
  printf("test_virial_lower_variance:\n");

  double Z = 2.0;
  double beta = 8.0;
  int P = RUNNING_ON_VALGRIND ? 32 : 256;
  int n_blocks = RUNNING_ON_VALGRIND ? 4 : 40;
  int samples_per_block = RUNNING_ON_VALGRIND ? 50 : 300;
  double tau = beta / P;
  int level = 4;

  pimc_result_t r =
      pimc_run(Z, P, tau, level, 400, n_blocks, samples_per_block, 777333ULL);

  printf("  error (thermodynamic)          = %.6f\n", r.error);
  printf("  error (virial)                 = %.6f\n", r.error_virial);
  printf("  variance ratio (thermo/virial) = %.2f\n",
         (r.error / r.error_virial) * (r.error / r.error_virial));

  check_true(r.error_virial < r.error,
             "virial estimator has lower standard error than the "
             "thermodynamic estimator on the same run");
}

static void test_pimc_run_invalid_input(void) {
  printf("test_pimc_run_invalid_input:\n");

  pimc_result_t r1 =
      pimc_run(2.0, 100, 0.05, 4, 10, 5, 10, 1ULL); // P not power of 2
  check_true(r1.n_blocks == 0, "non-power-of-2 P rejected");

  pimc_result_t r2 =
      pimc_run(2.0, 64, 0.05, 10, 10, 5, 10, 1ULL); // level > log2(P)
  check_true(r2.n_blocks == 0, "level > log2(P) rejected");
}

// NOTE: n_replicas=1 must reproduce pimc_run exactly, stream 0 has zero
// rng_jump() calls applied, and both route through same pimc_run_with_rng()
// sampling loop.
static void test_pimc_run_parallel_matches_serial_at_one_replica(void) {
  printf("test_pimc_run_parallel_matches_serial_at_one_replica:\n");

  double Z = 2.0;
  int P = RUNNING_ON_VALGRIND ? 16 : 64;
  double tau = 8.0 / P;
  int level = 3;

  pimc_result_t serial = pimc_run(Z, P, tau, level, 100, 10, 100, 606ULL);
  pimc_result_t parallel =
      pimc_run_parallel(1, Z, P, tau, level, 100, 10, 100, 606ULL);

  check_close(parallel.energy, serial.energy, 1e-12,
              "n_replicas=1 energy bit-identical to pimc_run");
  check_close(parallel.energy_virial, serial.energy_virial, 1e-12,
              "n_replicas=1 energy_virial bit-identical to pimc_run");
}

// Multiple independent PIMC ring-polymer walkers combined must still land close
// to exact helium ground state, with a sane nonzero inter-replica error bar.
// Lighter (P, n_blocks) than test_pimc_run_helium_accuracy to keep n_replicas *
// work bounded.
static void test_pimc_run_parallel_helium(void) {
  printf("test_pimc_run_parallel_helium:\n");

  int n_replicas = RUNNING_ON_VALGRIND ? 2 : 6;
  double Z = 2.0;
  int P = RUNNING_ON_VALGRIND ? 32 : 256;
  int n_blocks = RUNNING_ON_VALGRIND ? 4 : 20;
  int samples_per_block = RUNNING_ON_VALGRIND ? 50 : 200;
  double tau = 8.0 / P;
  int level = 4;

  pimc_result_t r = pimc_run_parallel(n_replicas, Z, P, tau, level, 400,
                                      n_blocks, samples_per_block, 24601ULL);

  printf("  parallel He (n_replicas=%d): E=%.6f +- %.6f  E_virial=%.6f +- %.6f "
         " n_blocks=%d\n",
         n_replicas, r.energy, r.error, r.energy_virial, r.error_virial,
         r.n_blocks);

  double E_exact = -2.9037;

  double tol = RUNNING_ON_VALGRIND ? 2.0 : 0.5;
  check_close(r.energy, E_exact, tol,
              "parallel PIMC energy lands close to exact helium ground state");
  check_true(r.error > 0.0, "parallel PIMC inter-replica error is nonzero");
  check_true(r.n_blocks == n_blocks * n_replicas,
             "n_blocks reports n_blocks_per_replica * n_replicas");
}

int main(void) {
  test_kelbg_fixtures();
  test_kelbg_limits();
  test_walker_alloc_and_init();
  test_bisection_invalid_input();
  test_pimc_run_invalid_input();
  test_pimc_run_helium_accuracy();
  test_virial_matches_thermodynamic();
  test_virial_lower_variance();
  test_pimc_run_parallel_matches_serial_at_one_replica();
  test_pimc_run_parallel_helium();

  if (failures == 0) {
    printf("\nAll test_pimc checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
