/*
 * Finite-system DMRG sweeps: Heisenberg antiferromagnetic spin chain
 *
 * physics/finite_dmrg.c extends infinite-system algorithm in dmrg.c
 * (Refernce: White, PRL 69, 2863 (1992)) with finite-size sweeping
 * procedure from White, PRB 48, 10345 (1993)
 * This example demonstrates its 2 defining properties, both direct consequences
 * of sweeping being able to rebuild every block using a better environment than
 * infinite algorithm ever had access to:
 *  1. For a fixed, truncation-limited bond dimension m, central-bond ground
 *     energy checkpoint improves (or holds) with every successive sweep
 *  2. At small m, the converged finite-sweep energy matches or improves on
 *     infinite algorithm's result for same (N, m)
 */

#include "../export/plot.h"
#include "../physics/dmrg.h"
#include "../physics/finite_dmrg.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Finite-system DMRG sweeps: Heisenberg antiferromagnetic "
         "chain\n\n");

  // NOTE: isotropic (critical) Heisenberg AFM - the case for a fixed, small
  // bond dimension since a gapless chain has no true "lossless" m
  double Jz = 1.0;
  double Jxy = 1.0;

  // Part 1: sweep-by-sweep convergence at fixed, small N and m
  int N_fixed = 20;
  int m_small = 8;
  int n_sweeps = 6;

  printf("  N=%d, m=%d fixed, tracking central-bond energy across sweeps:\n",
         N_fixed, m_small);
  finite_dmrg_result_t *r =
      finite_dmrg_run(N_fixed, Jz, Jxy, m_small, n_sweeps);
  if (!r) {
    fprintf(stderr, "finite_dmrg_run failed\n");

    return 1;
  }

  double *sweep_idx = malloc((size_t)n_sweeps * sizeof *sweep_idx);
  for (int s = 0; s < n_sweeps; s++) {
    sweep_idx[s] = s + 1;
    printf("    sweep %d:  E=%.10f\n", s + 1, r->sweep_energy[s]);
  }

  plot_opts_t opts1 = {0};
  opts1.title = "Finite DMRG central-bond energy vs sweep number (N=20, m=8)";
  opts1.xlabel = "sweep number";
  opts1.ylabel = "ground energy at N/2 bipartition";
  plot_line("finite_dmrg_sweep_convergence", PLOT_FORMAT_PNG, sweep_idx,
            r->sweep_energy, n_sweeps, &opts1);
  printf("\n   Saved finite_dmrg_sweep_convergence.png (should decrease "
         "monotonically then plateau)\n\n");

  free(sweep_idx);
  finite_dmrg_result_free(r);

  // Part 2: finite sweeps vs infinite algorithm at small, truncation-limited
  // m, across chain length
  int N_vals[] = {12, 16, 20, 24, 28};
  int n_N = (int)(sizeof(N_vals) / sizeof(N_vals[0]));
  int m_compare = 6;
  double *N_vals_d = malloc((size_t)n_N * sizeof *N_vals_d);
  double *e_infinite = malloc((size_t)n_N * sizeof *e_infinite);
  double *e_finite = malloc((size_t)n_N * sizeof *e_finite);

  printf("  m=%d (small, truncation-limited), sweeping chain length N:\n",
         m_compare);
  for (int i = 0; i < n_N; i++) {
    dmrg_result_t *r_inf = dmrg_run(N_vals[i], Jz, Jxy, m_compare);
    finite_dmrg_result_t *r_fin =
        finite_dmrg_run(N_vals[i], Jz, Jxy, m_compare, 5);
    if (!r_inf || !r_fin) {
      fprintf(stderr, "run failed at N=%d\n", N_vals[i]);
      free(r_inf);
      finite_dmrg_result_free(r_fin);

      continue;
    }

    N_vals_d[i] = N_vals[i];
    e_infinite[i] = r_inf->energy_per_site;
    e_finite[i] = r_fin->energy_per_site;

    printf("    N=%2d  E/N infinite=%.8f  E/N finite=%.8f  improvement=%.2e\n",
           N_vals[i], r_inf->energy_per_site, r_fin->energy_per_site,
           r_inf->energy_per_site - r_fin->energy_per_site);

    free(r_inf);
    finite_dmrg_result_free(r_fin);
  }

  plot_opts_t opts2 = {0};
  opts2.title = "Finite sweeps vs infinite algorithm at fixed small m=6";
  opts2.xlabel = "chain length N";
  opts2.ylabel = "E0/N";

  const double *series[2] = {e_infinite, e_finite};
  const char *labels[2] = {"infinite algorithm", "finite sweeps (5 sweeps)"};

  plot_lines("finite_dmrg_vs_infinite", PLOT_FORMAT_PNG, N_vals_d, series, 2,
             n_N, labels, &opts2);
  printf("\n   Saved finite_dmrg_vs_infinite.png (finite-sweep curve should "
         "sit at or below the infinite-algorithm curve: sweeping only ever "
         "re-optimizes previously-built blocks using more information)\n");

  free(N_vals_d);
  free(e_infinite);
  free(e_finite);

  return 0;
}
