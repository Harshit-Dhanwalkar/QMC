/*
 * Helium Ground State (Variational Method)
 *
 * Sweeps trial energy E(Z') for helium (and a couple of
 * helium-like ions)
 * variational minimum against the analytic Z'_opt = Z - 5/16.
 */

#include "../core/constants.h"
#include "../core/utils.h"
#include "../export/plot.h"
#include "../physics/helium.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void run_ion(double Z, const char *name) {
  double Zeff_opt;
  double E_hartree = helium_ground_state_energy_numeric(Z, 1e-12, &Zeff_opt);
  double E_eV = E_hartree * (AU_ENERGY / E_CHARGE);
  printf("   %-8s Z=%.0f:  Z'_opt=%.4f  E=%.6f Hartree (%.2f eV)\n", name, Z,
         Zeff_opt, E_hartree, E_eV);
}

int main(void) {
  printf(" > Helium Ground State (Effective-Charge Variational Method)\n\n");

  run_ion(1.0, "H-");
  run_ion(2.0, "He");
  run_ion(3.0, "Li+");
  run_ion(4.0, "Be2+");
  printf("\n");

  double Z = 2.0;
  int N = 200;
  double *Zeff = linspace(0.5, 3.0, N);
  double *E = malloc(N * sizeof *E);
  if (!Zeff || !E) {
    free(Zeff);
    free(E);
    return 1;
  }

  for (int i = 0; i < N; i++)
    E[i] = helium_variational_energy(Zeff[i], Z);

  double Zeff_analytic = helium_optimal_zeff_analytic(Z);
  double E_analytic = helium_ground_state_energy_analytic(Z);
  printf("   Helium: E(Z') minimized at Z'=%.4f, E=%.6f Hartree\n",
         Zeff_analytic, E_analytic);
  printf("   (theoritical: Z'=1.6875, E=-2.84765625 Hartree ~ -77.5 eV)\n\n");

  plot_opts_t opts = {0};
  opts.title = "Helium: E(Z') variational energy";
  opts.xlabel = "Z'";
  opts.ylabel = "E (Hartree)";
  // TODO: use save_wavefunction
  plot_line("helium_energy_curve", PLOT_FORMAT_PNG, Zeff, E, N, &opts);
  printf("   Saved helium_energy_curve.png (minimum near Z'=1.6875)\n");

  free(Zeff);
  free(E);

  return 0;
}
