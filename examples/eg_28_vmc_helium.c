/*
 * Variational Monte Carlo: Helium Ground State
 *
 * Sweeps Jastrow parameter b to show VMC energy landscape.
 * Finds optimal b via golden-section minimization.
 * Compare against plain product-orbital result (helium.c) and exact
 * ground-state energy.
 */

#include "../core/constants.h"
#include "../core/utils.h"
#include "../export/plot.h"
#include "../physics/helium.h"
#include "../physics/vmc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Variational Monte Carlo: Helium Ground State\n\n");

  double Zeff = 2.0;
  double E_exact = -2.9037;
  double E_product_orbital = helium_ground_state_energy_analytic(Zeff);

  printf("   Trial wavefunction: Psi_T = exp(-Z'(r1+r2)) * "
         "exp(r12/(2(1+b*r12)))\n");
  printf("   Z' fixed at bare nuclear charge (%.1f); sweeping Jastrow b:\n\n",
         Zeff);

  int N = 25;
  double b_min = 0.02, b_max = 0.5;
  double *b_vals = malloc(N * sizeof *b_vals);
  double *E_vals = malloc(N * sizeof *E_vals);
  if (!b_vals || !E_vals) {
    free(b_vals);
    free(E_vals);

    return 1;
  }

  for (int i = 0; i < N; i++) {
    double b = b_min + i * (b_max - b_min) / (N - 1);
    vmc_result_t r = vmc_run(Zeff, b, 500, 20000, 200, 0.9, 0.9, 1000ULL + i);
    b_vals[i] = b;
    E_vals[i] = r.mean;

    printf("   b=%.3f   E=%.5f +- %.5f Hartree\n", b, r.mean, r.error);
  }

  double b_opt;
  double E_opt = vmc_optimize_b(Zeff, b_min, b_max, 1000, 50000, 0.9, 0.9,
                                999ULL, 1e-3, &b_opt);

  printf("\n   Optimized: b_opt=%.4f   E=%.6f Hartree\n", b_opt, E_opt);
  printf("   Comparison:\n");
  printf("     Plain product orbital (no Jastrow): E=%.6f Hartree\n",
         E_product_orbital);
  printf("     VMC with optimized Jastrow:         E=%.6f Hartree\n", E_opt);
  printf("     Exact:                              E=%.6f Hartree\n", E_exact);
  printf("   Jastrow correlation recovers %.1f%% of the correlation energy "
         "this ansatz can reach.\n\n",
         100.0 * (E_product_orbital - E_opt) / (E_product_orbital - E_exact));

  plot_opts_t opts = {0};
  opts.title = "VMC Helium: E(b) Jastrow parameter sweep";
  opts.xlabel = "b";
  opts.ylabel = "E (Hartree)";
  plot_line("vmc_helium_energy_curve", PLOT_FORMAT_PNG, b_vals, E_vals, N,
            &opts);
  printf("   Saved vmc_helium_energy_curve.png (minimum near b=%.3f)\n", b_opt);

  free(b_vals);
  free(E_vals);

  return 0;
}
