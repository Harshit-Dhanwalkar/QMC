/*
 * Infinite-system DMRG: Heisenberg antiferromagnetic spin chain
 *
 * physics/dmrg.c implements White's infinite-system Density Matrix
 * Renormalization Group algorithm for open-boundary spin-1/2 XXZ chain.
 * This example demonstrates its two defining convergence properties:
 *  1. For fixed chain length N, ground-energy error (and discarded
 *     reduced-density-matrix weight, "truncation_error") shrinks rapidly as
 *     retained bond dimension m grows, becoming lossless once m is
 *     large enough to span chain's actual entanglement.
 *  2. For fixed (generous) m, per-site ground energy approaches
 *     exact Bethe-ansatz thermodynamic-limit value E0/N -> 1/4 - ln(2)
 *     (Reference: Hulthen 1938) as chain length N grows - isotropic Heisenberg
 *     chain is gapless/critical, so finite chains approach infinite-chain
 *     energy density only algebraically slowly in N.
 */

#include "../export/plot.h"
#include "../physics/dmrg.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Infinite-system DMRG: Heisenberg antiferromagnetic chain\n\n");

  // isotropic (critical) Heisenberg AFM
  double Jz = 1.0;
  double Jxy = 1.0;

  // Part 1: truncation-error convergence at fixed N
  int N_fixed = 16;
  int m_vals[] = {2, 4, 6, 8, 12, 16, 20, 24};
  int n_m = (int)(sizeof(m_vals) / sizeof(m_vals[0]));
  double *m_vals_d = malloc((size_t)n_m * sizeof *m_vals_d);
  double *trunc_err = malloc((size_t)n_m * sizeof *trunc_err);

  printf("  N=%d fixed, sweeping bond dimension m:\n", N_fixed);
  for (int i = 0; i < n_m; i++) {
    dmrg_result_t *r = dmrg_run(N_fixed, Jz, Jxy, m_vals[i]);
    if (!r) {
      fprintf(stderr, "dmrg_run failed at m=%d\n", m_vals[i]);

      continue;
    }

    m_vals_d[i] = m_vals[i];
    trunc_err[i] = r->truncation_error > 0 ? r->truncation_error : 1e-16;

    printf("    m=%2d  E=%.8f  E/N=%.8f  truncation_error=%.3e\n", m_vals[i],
           r->energy, r->energy_per_site, r->truncation_error);

    free(r);
  }

  plot_opts_t opts1 = {0};
  opts1.title = "DMRG truncation error vs bond dimension (N=16 XXZ chain)";
  opts1.xlabel = "bond dimension m";
  opts1.ylabel = "discarded density-matrix weight";
  plot_line("dmrg_truncation_error", PLOT_FORMAT_PNG, m_vals_d, trunc_err, n_m,
            &opts1);
  printf("\n   Saved dmrg_truncation_error.png (should fall steeply, "
         "flattening near machine precision once m spans entanglement)\n\n");

  free(m_vals_d);
  free(trunc_err);

  /* Part 2: approach to Bethe-ansatz thermodynamic limit */
  double bethe_limit =
      0.25 - log(2.0); // Reference: Hulthen 1938, exact for N->\infty
  int N_vals[] = {8, 12, 16, 20, 24, 30};
  int n_N = (int)(sizeof(N_vals) / sizeof(N_vals[0]));
  double *N_vals_d = malloc((size_t)n_N * sizeof *N_vals_d);
  double *e_per_site = malloc((size_t)n_N * sizeof *e_per_site);
  double *bethe_ref = malloc((size_t)n_N * sizeof *bethe_ref);
  int m_generous = 24;

  printf("  m=%d (generous), sweeping chain length N:\n", m_generous);
  for (int i = 0; i < n_N; i++) {
    dmrg_result_t *r = dmrg_run(N_vals[i], Jz, Jxy, m_generous);
    if (!r) {
      fprintf(stderr, "dmrg_run failed at N=%d\n", N_vals[i]);

      continue;
    }

    N_vals_d[i] = N_vals[i];
    e_per_site[i] = r->energy_per_site;
    bethe_ref[i] = bethe_limit;

    printf(
        "    N=%2d  E/N=%.8f  |E/N - (1/4-ln2)|=%.2e  truncation_error=%.2e\n",
        N_vals[i], r->energy_per_site, fabs(r->energy_per_site - bethe_limit),
        r->truncation_error);

    free(r);
  }

  printf("    exact thermodynamic limit 1/4 - ln(2) = %.8f\n\n", bethe_limit);

  plot_opts_t opts2 = {0};
  opts2.title = "DMRG per-site energy vs exact Bethe-ansatz limit";
  opts2.xlabel = "chain length N";
  opts2.ylabel = "E0/N";

  const double *series[2] = {e_per_site, bethe_ref};
  const char *labels[2] = {"DMRG (m=24)", "exact limit (1/4 - ln 2)"};

  plot_lines("dmrg_bethe_ansatz_limit", PLOT_FORMAT_PNG, N_vals_d, series, 2,
             n_N, labels, &opts2);

  printf("   Saved dmrg_bethe_ansatz_limit.png (DMRG curve should approach "
         "constant exact-limit line as N grows; convergence is algebraically "
         "slow since this chain is gapless/critical)\n");

  free(N_vals_d);
  free(e_per_site);
  free(bethe_ref);

  return 0;
}
