/*
 * Path Integral Monte Carlo: Helium Ground State
 *
 * Finite-temperature PIMC (Kelbg-regularized Coulomb, bisection sampling) at
 * low temperature (large \beta) projects toward the ground-state energy, same
 * target as VMC/DMC but via a completely different method: no trial
 * wavefunction, direct sampling of the imaginary-time path integral.
 */

#include "../core/constants.h"
#include "../export/plot.h"
#include "../physics/pimc.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Path Integral Monte Carlo: Helium Ground State\n\n");

  double Z = 2.0;
  double beta = 8.0; // low enough temperature that the finite-T PIMC estimator
                     // approaches the ground state energy
  double E_exact = -2.9037;

  printf("   Two electrons as distinguishable ring polymers (valid here: "
         "helium ground state is nodeless, so there's no fermion sign problem "
         "for this system specifically).\n");
  printf("   Kelbg-regularized Coulomb pair potentials (finite at r=0) avoid "
         "the path-collapse catastrophe a bare Coulomb singularity causes; "
         "bisection moves resample whole path segments from the exact "
         "free-particle distribution instead of diffusing bead by bead.\n\n");

  printf("   \\beta=%.1f fixed; \\tau-convergence via P: 64 -> 512 (\\tau = "
         "\\beta / P):\n\n",
         beta);

  const int Ps[] = {64, 128, 256, 512};
  int n = 4;
  double *tau_vals = malloc(n * sizeof *tau_vals);
  double *E_vals = malloc(n * sizeof *E_vals);
  double *E_virial_vals = malloc(n * sizeof *E_virial_vals);
  if (!tau_vals || !E_vals || !E_virial_vals) {
    free(tau_vals);
    free(E_vals);
    free(E_virial_vals);

    return 1;
  }

  printf(
      "   Two estimators are computed on the same sampled configurations:\n");
  printf("  the thermodynamic estimator (E) and the lower-variance virial "
         "estimator (E_virial, same expectation value, derived via a "
         "coordinate-rescaling identity\n");

  for (int i = 0; i < n; i++) {
    int P = Ps[i];
    double tau = beta / P;
    int level = 4; // segment length 16 beads

    pimc_result_t r = pimc_run(Z, P, tau, level, 400, 30, 250, 20260802ULL);

    tau_vals[i] = tau;
    E_vals[i] = r.energy;
    E_virial_vals[i] = r.energy_virial;

    printf("   P=%4d  \\tau=%.5f  E=%.6f+-%.6f   E_virial=%.6f+-%.6f  "
           "acceptance=%.3f\n",
           P, tau, r.energy, r.error, r.energy_virial, r.error_virial,
           r.acceptance_rate);
  }

  printf("\n   Exact ground state energy: %.6f Hartree\n", E_exact);
  printf("   PIMC virial estimator (P=512) is %.4f Hartree from exact: "
         "consistent\n",
         E_virial_vals[n - 1] - E_exact);

  plot_opts_t opts = {0};
  opts.title = "PIMC Helium: E(\\tau) convergence at fixed \\beta=8";
  opts.xlabel = "\\tau";
  opts.ylabel = "E_virial (Hartree)";
  plot_line("pimc_helium_tau_convergence", PLOT_FORMAT_PNG, tau_vals,
            E_virial_vals, n, &opts);
  printf("   Saved pimc_helium_tau_convergence.png\n");

  free(tau_vals);
  free(E_vals);
  free(E_virial_vals);

  return 0;
}
