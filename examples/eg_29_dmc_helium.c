/*
 * Diffusion Monte Carlo: Helium Ground State
 *
 * DMC projecting past VMC's variational bound toward exact ground state energy.
 * For system trial wavefunction has no nodes, so plain DMC reaches exact answer
 * up to finite-timestep and finite-population effects (no fixed-node
 * approximation needed.)
 */

#include "../core/constants.h"
#include "../physics/dmc.h"
#include "../physics/helium.h"
#include "../physics/vmc.h"
#include <stdio.h>

int main(void) {
  printf(" > Diffusion Monte Carlo: Helium Ground State\n\n");

  double Zeff = 2.0;
  double b = 0.15;
  double E_exact = -2.9037;
  double E_product_orbital = helium_ground_state_energy_analytic(Zeff);

  printf("   Trial wavefunction: Slater-Jastrow (Z'=%.1f, b=%.2f)\n\n", Zeff,
         b);

  printf("   Running VMC (200000 samples)...\n");
  vmc_result_t vmc_r = vmc_run(Zeff, b, 2000, 200000, 200, 0.9, 0.9, 1234ULL);
  printf("     E = %.6f +- %.6f Hartree\n\n", vmc_r.mean, vmc_r.error);

  printf("   Running DMC (population ~500, \\tau=0.01, 30 blocks x 200 "
         "generations)...\n");
  dmc_result_t dmc_r =
      dmc_run(Zeff, b, 500, 1500, 0.01, 1000, 30, 200, 5678ULL);
  printf("     mixed estimator:  E = %.6f +- %.6f Hartree\n",
         dmc_r.energy_mixed, dmc_r.error_mixed);
  printf("     growth estimator: E = %.6f +- %.6f Hartree\n",
         dmc_r.energy_growth, dmc_r.error_growth);
  printf("     mean population=%.1f   acceptance=%.3f\n\n",
         dmc_r.mean_population, dmc_r.acceptance_rate);

  printf("   Comparison (Hartree):\n");
  printf("     Plain product orbital (no correlation): %.6f\n",
         E_product_orbital);
  printf("     VMC (Jastrow-correlated, this run):     %.6f\n", vmc_r.mean);
  printf("     DMC (this run):                         %.6f\n",
         dmc_r.energy_mixed);
  printf("     Exact:                                  %.6f\n\n", E_exact);

  double vmc_gap = vmc_r.mean - E_exact;
  double dmc_gap = dmc_r.energy_mixed - E_exact;
  printf("   VMC is %.4f Hartree above exact; DMC is %.4f Hartree above "
         "exact.\n",
         vmc_gap, dmc_gap);
  printf("   DMC recovers the correlation energy VMC's finite Jastrow "
         "ansatz misses.\n");

  return 0;
}
