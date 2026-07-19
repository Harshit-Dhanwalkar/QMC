/*
 * Two-Level Rabi Oscillations
 *
 * Plots the excited-state probability P_e(t) for a driven two-level
 * system, on resonance and detuned, and verifies the resonant
 * pi-pulse gives full population inversion.
 */

#include "../core/complex.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../export/plot.h"
#include "../physics/rabi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Two-Level Rabi Oscillations\n\n");

  double Omega = 1.0;
  int N = 300;
  double t_max = 4.0 * M_PI / Omega; // a few full resonant periods
  double *t = linspace(0.0, t_max, N);
  if (!t)
    return 1;

  double *P_resonant = malloc(N * sizeof *P_resonant);
  double *P_detuned = malloc(N * sizeof *P_detuned);
  if (!P_resonant || !P_detuned) {
    free(t);
    free(P_resonant);
    free(P_detuned);
    return 1;
  }

  double Delta_detuned = 2.0 * Omega;
  for (int i = 0; i < N; i++) {
    P_resonant[i] = rabi_excited_probability(t[i], Omega, 0.0);
    P_detuned[i] = rabi_excited_probability(t[i], Omega, Delta_detuned);
  }

  printf("   \\Omega=%.2f, resonant Delta=0, detuned Delta=%.2f\n\n", Omega,
         Delta_detuned);

  double t_pi = M_PI / Omega;
  printf("   Resonant pi-pulse at t=%.4f: P_e=%.6f (expect 1.0)\n\n", t_pi,
         rabi_excited_probability(t_pi, Omega, 0.0));

  const double *ys[2] = {P_resonant, P_detuned};
  const char *labels[2] = {"resonant", "detuned"};
  plot_opts_t opts = {0};

  opts.title = "Rabi Oscillations: P_excited(t)";
  opts.xlabel = "t";
  opts.ylabel = "P_e";

  plot_lines("rabi_oscillations", PLOT_FORMAT_PNG, t, ys, 2, N, labels, &opts);
  printf("   Saved rabi_oscillations.png\n");

  free(t);
  free(P_resonant);
  free(P_detuned);
  return 0;
}
