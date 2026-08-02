/*
Test: Path Integral Monte Carlo for helium ground state (Kelbg-regularized
Coulomb, bisection sampling).

1. kelbg_potential / kelbg_energy_correction at fixed (r, \lambda, q) must match
independently-derived closed-form values - deterministic, no MC noise.
2. Kelbg potential must be finite at r=0 and must reduce to the bare Coulomb
potential q/r as \lambda -> 0 (\tau -> 0 high-temperature limit).
3. Zero-charge sanity check: q=0 must give exactly V=0 for both functions, at
any r and \lambda (confirms charge scaling isn't accidentally offset).
4. Invalid-input / edge-case handling for allocation and bisection moves.
5. Full pimc_run output close to exact helium ground state energy (-2.9037
Hartree) at validated (P, \tau, level) parameters, with reasonable acceptance
rate.
*/

// FIX: Thermodynamic energy estimator's variance grows with P (PIMC issue:
// kinetic term is difference of two large, partially-cancelling quantities), so
// pushing P much beyond ~512 at fixed statistics needs proportionally more
// samples for tight comparison, or a virial estimator (which is not
// implemented) to avoid growing variance altogether.
// HACK: This test uses P=512 specifically because it's validated and found
// during development, not largest P that runs.

#include "../core/random.h"
#include "../physics/pimc.h"
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

static void test_kelbg_fixtures(void) {
  printf("test_kelbg_fixtures:\n");

  {
    double r = 0.618, lambda = sqrt(0.139 / 0.5), q = 1.0; // tau=0.139
    double v = kelbg_potential(r, lambda, q);
    double vc = kelbg_energy_correction(r, lambda, q);
    check_close(v, 1.5359326977018246, 1e-9, "kelbg_potential fixture 1");
    check_close(vc, 1.3722250340689774, 1e-9,
                "kelbg_energy_correction fixture 1");
  }
  {
    double r = 2.885, lambda = sqrt(0.365 / 0.5), q = -3.0; // tau=0.365
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

// Full run at validated (P, \\tau, level) parameters
static void test_pimc_run_helium_accuracy(void) {
  printf("test_pimc_run_helium_accuracy:\n");

  double Z = 2.0;
  double beta = 8.0;
  int P = 512;
  double tau = beta / P;
  int level = 4;

  pimc_result_t r = pimc_run(Z, P, tau, level, 400, 40, 300, 909090ULL);

  printf("  P=%d tau=%.5f beta=%.2f: E=%.6f +- %.6f  acceptance=%.4f  "
         "n_blocks=%d\n",
         P, tau, P * tau, r.energy, r.error, r.acceptance_rate, r.n_blocks);

  double E_exact = -2.9037;

  check_true(r.n_blocks == 40, "requested number of blocks completed");
  check_true(r.acceptance_rate > 0.3,
             "bisection acceptance stays well above the near-zero "
             "pathology a naive full-ring move shows at this P");

  check_close(r.energy, E_exact, 0.4,
              "PIMC energy lands close to exact helium ground state");
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

int main(void) {
  test_kelbg_fixtures();
  test_kelbg_limits();
  test_walker_alloc_and_init();
  test_bisection_invalid_input();
  test_pimc_run_invalid_input();
  test_pimc_run_helium_accuracy();

  if (failures == 0) {
    printf("\nAll test_pimc checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
