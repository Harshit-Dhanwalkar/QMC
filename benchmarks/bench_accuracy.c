/*
 * Benchmark: VMC / DMC ground-state energy vs. published nonrelativistic
 * two-electron reference values.
 *
 * physics/vmc.c and physics/dmc.c implement Slater-Jastrow VMC/DMC for two-electron atoms/ions only (He, H-, Li+, Be2+, ...).
 *
 * Reference: Published values are clamped-nucleus, nonrelativistic ground-state energies (Hartree) from the Pekeris-type high-precision literature for the helium isoelectronic sequence.
 *
 * USAGE:
 *   make PLOT_BACKEND=NONE SANITIZE=0 build/bench_accuracy
 *   ./build/bench_accuracy
 */

#include "../physics/dmc.h"
#include "../physics/vmc.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>

int main(void) {
  printf("# QMC Accuracy vs. Published Results (He isoelectronic sequence, "
         "two-electron only)\n");
  printf("# System | VMC Energy | DMC Energy | Published | VMC Error | DMC "
         "Error\n");

  struct {
    char name[16];
    double Z;
    double Zeff_guess; /* starting point for b-optimization; Z minus the
                        * standard variational electron-screening constant
                        * 5/16 */
    double published;  /* Hartree, clamped-nucleus nonrelativistic */
  } systems[] = {
      {"He", 2.0, 2.0 - 5.0 / 16.0, -2.903724},
      {"Li+", 3.0, 3.0 - 5.0 / 16.0, -7.279913},
      {"Be2+", 4.0, 4.0 - 5.0 / 16.0, -13.655566},
  };
  int n_systems = 3;

  for (int i = 0; i < n_systems; i++) {
    double Z = systems[i].Z;
    double Zeff = systems[i].Zeff_guess;

    /* Step sizes and DMC timestep shrink with Z: electrons localize closer
     * to a more charged nucleus (~1/Z length scale), so a fixed step size
     * tuned for He would over-propose (low acceptance) for Li+/Be2+. */
    double step = 0.9 / Zeff;
    double tau = 0.01 / (Zeff * Zeff);

    double b_opt;
    double e_vmc_opt =
        vmc_optimize_b(Z, Zeff, 0.02, 0.6, 500, 20000, step, step,
                       1000ULL + (uint64_t)i, 1e-3, &b_opt);

    /* Re-run VMC at the optimized b with more samples for a tighter
     * reported error bar than the optimizer's own inner evaluations use. */
    vmc_result_t vmc_r = vmc_run(Z, Zeff, b_opt, 2000, 100000, 200, step, step,
                                 2000ULL + (uint64_t)i);

    dmc_result_t dmc_r =
        dmc_run(Z, Zeff, b_opt, /*target_population=*/500,
                /*max_population=*/1500, tau, /*n_equilibration=*/1000,
                /*n_blocks=*/30, /*block_size=*/200, 3000ULL + (uint64_t)i);

    double vmc_err = fabs(vmc_r.mean - systems[i].published);
    double dmc_err = fabs(dmc_r.energy_mixed - systems[i].published);

    printf("%-6s %.6f %.6f %.6f %.6f %.6f   (b_opt=%.4f, VMC E=%.6f from "
           "optimizer)\n",
           systems[i].name, vmc_r.mean, dmc_r.energy_mixed,
           systems[i].published, vmc_err, dmc_err, b_opt, e_vmc_opt);
    fflush(stdout);
  }

  return 0;
}
