/*
Test: Diffusion Monte Carlo for helium ground state.

1. dmc_drift_velocity at fixed configurations (deterministic, no MC noise).
2. Full dmc_run must land close to exact helium ground state energy (-2.9037
   Hartree) - should be much closer than VMC's variational estimate (-2.84765625
   from plain product orbital, or ~-2.878 from VMC with same Jastrow).
3. Acceptance rate and mean population must be in a sane range, confirming
   sampler isn't degenerate and population control is holding near target.
*/

#include "../core/random.h"
#include "../physics/dmc.h"
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

static void test_drift_velocity_fixtures(void) {
  printf("test_drift_velocity_fixtures:\n");

  {
    vmc_walker_t w = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
    double b = 0.2;
    double d0[3], d1[3];

    dmc_drift_velocity(&w, 0, 2.0, b, d0);
    dmc_drift_velocity(&w, 1, 2.0, b, d1);
    check_close(d0[0], -1.7448979591836733, 1e-9, "drift0[0] case 1");
    check_close(d0[1], 0.0, 1e-9, "drift0[1] case 1");
    check_close(d1[0], 1.7448979591836733, 1e-9, "drift1[0] case 1");
  }
  {
    vmc_walker_t w = {{0.5, 0.3, -0.2}, {-0.4, 0.6, 0.1}};
    double b = 0.15;
    double d0[3], d1[3];

    dmc_drift_velocity(&w, 0, 2.0, b, d0);
    dmc_drift_velocity(&w, 1, 2.0, b, d1);
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
  dmc_drift_velocity(&w, 0, 2.0, 0.2, d);
  check_true(isfinite(d[0]) && isfinite(d[1]) && isfinite(d[2]) &&
                 d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.0,
             "r1->0 configuration returns finite zero-vector guard");
}

static void test_dmc_run_accuracy(void) {
  printf("test_dmc_run_accuracy:\n");

  double Z = 2.0, Zeff = 2.0, b = 0.15;
  dmc_result_t r = dmc_run(Z, Zeff, b,
                           /*target_population=*/200,
                           /*max_population=*/600,
                           /*\tau=*/0.01,
                           /*n_equilibration=*/500,
                           /*n_blocks=*/20,
                           /*block_size=*/200,
                           /*seed=*/20260731ULL);

  printf("  mixed:  %.6f +- %.6f Hartree\n", r.energy_mixed, r.error_mixed);
  printf("  growth: %.6f +- %.6f Hartree\n", r.energy_growth, r.error_growth);
  printf("  mean_population=%.2f  acceptance=%.4f  n_blocks=%d\n",
         r.mean_population, r.acceptance_rate, r.n_blocks);

  double E_exact = -2.9037;
  double E_vmc_jastrow = -2.878; // \aprrox value from VMC session's b-scan
  double E_product_orbital = helium_ground_state_energy_analytic(2.0);

  /* DMC should land close to exact (much tighter than a VMC-level tolerance
   * would be, but still allowing for residual timestep/population-control bias
   * at finite \tau and finite population). */
  check_close(r.energy_mixed, E_exact, 0.02,
              "mixed estimator close to exact ground state");

  check_true(r.energy_mixed < E_vmc_jastrow,
             "DMC energy improves on VMC-with-Jastrow estimate");
  check_true(r.energy_mixed < E_product_orbital,
             "DMC energy improves on plain product-orbital estimate");

  /* Shouldn't wildly undershoot exact value either */
  check_true(r.energy_mixed > E_exact - 0.05,
             "DMC energy does not undershoot exact value implausibly");

  check_true(r.acceptance_rate > 0.7 && r.acceptance_rate <= 1.0,
             "acceptance rate in reasonable range (expect high, ~0.9+, for "
             "Metropolis-corrected drift-diffusion at this \\tau)");

  check_true(fabs(r.mean_population - 200.0) < 20.0,
             "population control holds near target (200 +/- 20)");
}

// DMC for a different two-electron ion (Li+, Z=3): must land closer to exact
// reference than VMC's Li+ estimate does and must respect variational theorem.
static void test_dmc_run_different_ion(void) {
  printf("test_dmc_run_different_ion:\n");

  double Z = 3.0, Zeff = 2.6, b = 0.12; // Li+
  double E_exact_liplus = -7.2799133;

  vmc_result_t vmc_r =
      vmc_run(Z, Zeff, b, 1000, 100000, 200, 0.9, 0.9, 55221ULL);
  dmc_result_t dmc_r = dmc_run(Z, Zeff, b,
                               /*target_population=*/200,
                               /*max_population=*/600,
                               /*tau=*/0.01,
                               /*n_equilibration=*/500,
                               /*n_blocks=*/20,
                               /*block_size=*/200,
                               /*seed=*/55222ULL);

  printf("  Li+ VMC: %.6f +- %.6f Hartree\n", vmc_r.mean, vmc_r.error);
  printf("  Li+ DMC (mixed): %.6f +- %.6f Hartree\n", dmc_r.energy_mixed,
         dmc_r.error_mixed);

  check_true(dmc_r.energy_mixed >= E_exact_liplus - 3.0 * dmc_r.error_mixed,
             "Li+ DMC respects the variational theorem");
  check_true(dmc_r.energy_mixed < vmc_r.mean,
             "Li+ DMC improves on Li+ VMC estimate (projects toward exact)");
}

int main(void) {
  test_drift_velocity_fixtures();
  test_degenerate_guard();
  test_dmc_run_accuracy();
  test_dmc_run_different_ion();

  if (failures == 0) {
    printf("\nAll test_dmc checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
