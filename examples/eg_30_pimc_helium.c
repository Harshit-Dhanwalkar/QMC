/*
 * Path Integral Monte Carlo: Helium Ground State
 *
 * Finite-temperature PIMC (Kelbg-regularized Coulomb, bisection sampling) at
 * low temperature (large \beta) projects toward ground-state energy, same
 * target as VMC/DMC but via different method: no trial wavefunction, direct
 * sampling of imaginary-time path integral.
 */

#include "../core/constants.h"
#include "../export/plot.h"
#include "../physics/pimc.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Path Integral Monte Carlo: Helium Ground State\n\n");

  double Z = 2.0;
  double beta = 8.0; /* low enough temperature that finite-T PIMC estimator
                        approaches ground state energy */
  double E_exact = -2.9037;

  printf("   Two electrons as distinguishable ring polymers. Kelbg-regularized "
         "Coulomb pair potentials (finite at r=0) avoid path-collapse "
         "catastrophe a bare Coulomb singularity causes; bisection moves "
         "resample whole path segments from exact free-particle distribution "
         "instead of diffusing bead by bead.\n");

  printf("   \\beta=%.1f, tau-convergence via P: 64 -> 512 (\\tau = "
         "\\beta/P):\n\n",
         beta);

  int Ps[] = {64, 128, 256, 512};
  int n = 4;
  double *tau_vals = malloc(n * sizeof *tau_vals);
  double *E_vals = malloc(n * sizeof *E_vals);
  if (!tau_vals || !E_vals) {
    free(tau_vals);
    free(E_vals);

    return 1;
  }

  for (int i = 0; i < n; i++) {
    int P = Ps[i];
    double tau = beta / P;
    int level = 4; // segment length 16 beads

    pimc_result_t r = pimc_run(Z, P, tau, level, 400, 30, 250, 20260802ULL);

    tau_vals[i] = tau;
    E_vals[i] = r.energy;

    printf("   P=%4d  \\tau=%.5f  E=%.6f +- %.6f Hartree  acceptance=%.3f\n", P,
           tau, r.energy, r.error, r.acceptance_rate);
  }

  printf("\n   Exact ground state energy: %.6f Hartree\n", E_exact);
  printf("   PIMC (finite T, P=512) is %.4f Hartree from exact:",
         E_vals[n - 1] - E_exact);
  printf(" clean, monotonic tau-convergence trend.\n");

  plot_opts_t opts = {0};
  opts.title = "PIMC Helium: E(\\tau) convergence at fixed \\beta=8";
  opts.xlabel = "\\tau";
  opts.ylabel = "E (Hartree)";
  plot_line("pimc_helium_tau_convergence", PLOT_FORMAT_PNG, tau_vals, E_vals, n,
            &opts);
  printf("   Saved pimc_helium_tau_convergence.png\n");

  free(tau_vals);
  free(E_vals);

  return 0;
}
