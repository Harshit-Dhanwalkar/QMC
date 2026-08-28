/*
 * Benchmark: VMC statistical convergence vs. sample count, He atom.
 *
 * Real entry point is vmc_run(Z, Zeff, b, n_equilibration, n_samples,
 * block_size, step1, step2, seed) - see physics/vmc.h. "n_samples" there is
 * sweeps per single chain (this file's "walkers" axis), not a walker-count
 * parameter
 * VMC here is single-walker-chain Metropolis, not a population method
 * (that's DMC's dmc_population_t).
 *
 * b is fixed at a reasonable pre-optimized value (see eg_28_vmc_helium.c)
 * rather than re-optimized per sample count, so scan below isolates pure
 * statistical (1 / \sqrt(N)) convergence and doesn't conflate it with optimizer
 * noise at small sample counts.
 *
 * USAGE:
 *   make PLOT_BACKEND=NONE SANITIZE=0 build/bench_vmc_convergence
 *   ./build/bench_vmc_convergence
 */

#include "../physics/vmc.h"
#include <bits/time.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

static double wall_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);

  double time = ts.tv_sec + ts.tv_nsec * 1e-9;

  // return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
  return time;
}

int main(void) {
  printf("# VMC Convergence Benchmark: He Atom\n");
  printf("# Samples | Mean Energy | Reported Std Error | Wall Time (s) | "
         "Error vs. Exact\n");

  double Z = 2.0, Zeff = 2.0, b = 0.15;
  double exact_energy = -2.903724;
  int n_equilibration = 1000;

  const int sample_counts[] = {100, 500, 1000, 5000, 10000, 50000};
  int n_points = 6;

  for (int i = 0; i < n_points; i++) {
    int n_samples = sample_counts[i];
    /* At least ~5 blocks even at the smallest sample count, so the reported
     * block-averaged std error isn't degenerate (a single block reports 0). */
    int block_size = n_samples / 5 < 200 ? (n_samples / 5) : 200;
    if (block_size < 1) {
      block_size = 1;
    }

    double t0 = wall_seconds();
    vmc_result_t r = vmc_run(Z, Zeff, b, n_equilibration, n_samples, block_size,
                             0.9, 0.9, 1000ULL + (uint64_t)i);
    double elapsed = wall_seconds() - t0;

    double error = fabs(r.mean - exact_energy);

    printf("%d %f %f %f %f\n", n_samples, r.mean, r.error, elapsed, error);

    fflush(stdout);
  }

  printf("\n# Note: reported std error above is the block-averaged standard "
         "error on the mean (see vmc_result_t.error), which should shrink "
         "roughly as 1/sqrt(n_samples) if blocks are long enough to average "
         "out Metropolis autocorrelation : a cross-check against empirical "
         "|error vs exact| column.\n");

  return 0;
}
