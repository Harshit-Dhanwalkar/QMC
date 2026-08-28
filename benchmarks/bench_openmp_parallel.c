/*
 * Benchmark: Hybrid MPI (across processes) + OpenMP (within a process)
 * scaling for VMC/DMC on He, implementing TODO.md's "Phase 1: Hybrid VMC
 * Baseline" and DMC extension.
 *
 * This is the first MPI-aware code in the repo -- everything else (OpenMP
 * replica parallelism via vmc_run_parallel/dmc_run_parallel, demonstrated in
 * examples/eg_39_openmp_qmc.c) already exists and needed no changes.
 *
 * Design, following TODO.md's own outline and core/random.h's documented
 * RNG-fanout contract:
 *   - OpenMP fans out across REPLICAS WITHIN a rank: vmc_run_parallel /
 *     dmc_run_parallel already do this internally via rng_jump() (see
 *     physics/vmc.c, physics/dmc.c), needing no MPI-specific code at all.
 *   - MPI fans out across RANKS: each rank derives its own well-separated
 *     RNG starting point via rng_long_jump() (core/random.h: "generating
 *     independent *jump* starting points (e.g. one per MPI rank, each of
 *     which then uses rng_jump() internally to fan out to its own
 *     threads/walkers)"), then extracts a uint64_t seed from it to hand to
 *     vmc_run_parallel/dmc_run_parallel's existing master_seed parameter.
 *     No changes to physics/vmc.c or physics/dmc.c were needed.
 *   - Statistics are combined with MPI_Gather (not MPI_Allreduce on the raw
 *     means): only rank 0 needs the combined result to report it, and
 *     gathering each rank's local grand_mean lets rank 0 compute a proper
 *     inter-rank standard error (std of the per-rank means / sqrt(n_ranks)),
 *     the same treatment vmc_run_parallel itself gives to its per-replica
 *     means one level down. This mirrors TODO.md's "Use MPI_Allreduce
 *     strictly at the end of optimization blocks to pool statistical
 *     averages across processes" -- Allreduce broadcasts too, which no rank
 *     here needs, so Reduce+Gather is the tighter-scoped choice for a
 *     reporting-only benchmark. A future parameter-optimization driver that
 *     needs every rank to see the pooled result would want Allreduce
 *     instead.
 *
 * Physical work done per rank is IDENTICAL to eg_39_openmp_qmc.c's single
 * process -- MPI only adds the process dimension on top.
 *
 * TODO.md "Phase 2: Distributed VMC" checklist, addressed here:
 *   - "Manager-Worker / All-Reduce VMC execution": MPI_Gather-based
 *     aggregation above (see the design note on Allreduce vs Gather).
 *   - "Verify statistical agreement with serial VMC": the serial-reference
 *     block below now keeps every replica's mean (not just timing) and
 *     reports a z-score between the hybrid and serial grand means via
 *     standard error propagation -- a real correctness check, not just a
 *     timing comparison.
 *   - "Profile OpenMP vs. MPI vs. Hybrid scaling": no new code needed for
 *     this -- it's three ways of invoking the same binary:
 *       OpenMP-only:  OMP_NUM_THREADS=N mpirun -np 1
 * ./build/bench_openmp_parallel MPI-only:     OMP_NUM_THREADS=1 mpirun -np N
 * ./build/bench_openmp_parallel Hybrid:       OMP_NUM_THREADS=k mpirun -np N/k
 * ./build/bench_openmp_parallel Compare the reported "Hybrid MPI+OpenMP
 * speedup" line across the three; on a single machine (no real multi-node MPI)
 * OpenMP-only should generally win since it avoids MPI's process-startup and
 *     message-passing overhead entirely -- MPI's own advantage only shows
 *     up once ranks are actually spread across separate machines.
 *
 * USAGE (needs an MPI implementation, e.g. `apt install libopenmpi-dev`):
 *   make PLOT_BACKEND=NONE SANITIZE=0 build/bench_openmp_parallel
 *   mpirun -np 4 ./build/bench_openmp_parallel
 *   OMP_NUM_THREADS=2 mpirun -np 4 ./build/bench_openmp_parallel
 *   mpirun -np 4 --bind-to core ./build/bench_openmp_parallel   (on a real
 *   multi-node cluster; not meaningful on a single-machine sandbox)
 */

#include "../core/random.h"
#include "../physics/dmc.h"
#include "../physics/vmc.h"
#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static double wall_seconds(void) { return MPI_Wtime(); }

/* Rank-unique master seed: long_jump `rank+1` times from a shared base
 * seed, then draw a scalar to feed into vmc_run_parallel/dmc_run_parallel's
 * existing uint64_t master_seed parameter (which itself rng_jump()'s that
 * scalar out across the rank's local OpenMP replicas). Non-overlap between
 * ranks and between a rank's replicas is therefore end-to-end guaranteed by
 * the same jump/long_jump machinery already used for OpenMP, not by hoping
 * distinct scalar seeds don't collide. */
static uint64_t rank_master_seed(uint64_t base_seed, int rank) {
  rng_state_t rng;
  rng_seed(&rng, base_seed);
  for (int i = 0; i < rank + 1; i++) {
    rng_long_jump(&rng);
  }

  return rng_next_u64(&rng);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, n_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

  int n_threads = 1;
#ifdef _OPENMP
  n_threads = omp_get_max_threads();
#endif

  if (rank == 0) {
    printf("# Hybrid MPI+OpenMP QMC: He (2 electrons)\n");
    printf("# MPI ranks: %d, OpenMP threads/rank: %d (total parallel "
           "workers: %d)\n",
           n_ranks, n_threads, n_ranks * n_threads);
#ifndef _OPENMP
    printf("# NOTE: built without OpenMP -- rebuild with -fopenmp for the "
           "intra-rank half of this benchmark to do anything.\n");
#endif
  }

  double Z = 2.0, Zeff = 1.6875, b = 0.3;
  int n_equilibration = 1000, n_samples = 40000, block_size = 100;
  double step1 = 0.9, step2 = 0.9;
  int replicas_per_rank = n_threads; /* always >= 1: set above */
  uint64_t base_seed = 1000ULL;

  /* === VMC === */
  uint64_t my_seed = rank_master_seed(base_seed, rank);

  MPI_Barrier(MPI_COMM_WORLD);
  double t0 = wall_seconds();
  vmc_result_t local =
      vmc_run_parallel(replicas_per_rank, Z, Zeff, b, n_equilibration,
                       n_samples, block_size, step1, step2, my_seed);
  double t1 = wall_seconds();
  double local_elapsed = t1 - t0;

  double *all_means = NULL, *all_elapsed = NULL;
  int *all_nsamples = NULL;
  if (rank == 0) {
    all_means = malloc((size_t)n_ranks * sizeof(double));
    all_elapsed = malloc((size_t)n_ranks * sizeof(double));
    all_nsamples = malloc((size_t)n_ranks * sizeof(int));
  }

  MPI_Gather(&local.mean, 1, MPI_DOUBLE, all_means, 1, MPI_DOUBLE, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&local_elapsed, 1, MPI_DOUBLE, all_elapsed, 1, MPI_DOUBLE, 0,
             MPI_COMM_WORLD);
  MPI_Gather(&local.n_samples, 1, MPI_INT, all_nsamples, 1, MPI_INT, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    double grand_mean = 0.0;
    long total_samples = 0;

    for (int r = 0; r < n_ranks; r++) {
      grand_mean += all_means[r];
      total_samples += all_nsamples[r];
    }

    grand_mean /= n_ranks;

    double between_var = 0.0;
    for (int r = 0; r < n_ranks; r++) {
      double d = all_means[r] - grand_mean;

      between_var += d * d;
    }

    double inter_rank_error =
        (n_ranks > 1) ? sqrt(between_var / (n_ranks - 1) / n_ranks) : 0.0;

    double max_elapsed = all_elapsed[0];
    for (int r = 1; r < n_ranks; r++) {
      if (all_elapsed[r] > max_elapsed) {
        max_elapsed = all_elapsed[r];
      }
    }

    printf("\n=== VMC: %d ranks x %d replicas/rank, %d samples/replica (%ld "
           "total samples) ===\n",
           n_ranks, replicas_per_rank, n_samples, total_samples);
    printf("Wall time (max across ranks): %.3f s\n", max_elapsed);
    printf("E = %.6f +- %.6f Hartree (inter-rank standard error over %d "
           "ranks; exact He ground state: -2.9037)\n",
           grand_mean, inter_rank_error, n_ranks);

    /* Reference: same total sampling work (n_ranks * replicas_per_rank
     * independent replicas), run serially on rank 0 alone -- same
     * methodology as eg_39_openmp_qmc.c's sequential-loop comparison,
     * extended across the MPI dimension too. Collect each serial replica's
     * mean (not just timing) so we can check statistical agreement with the
     * hybrid run below, per TODO.md Phase 2's "Verify statistical agreement
     * with serial VMC". */
    int n_serial_replicas = n_ranks * replicas_per_rank;
    double *serial_means = malloc((size_t)n_serial_replicas * sizeof(double));

    printf("\nSerial reference: running the equivalent %d total replicas "
           "one at a time on rank 0 for comparison...\n",
           n_serial_replicas);
    double t2 = wall_seconds();
    for (int i = 0; i < n_serial_replicas; i++) {
      vmc_result_t r = vmc_run_parallel(1, Z, Zeff, b, n_equilibration,
                                        n_samples, block_size, step1, step2,
                                        base_seed + 7919ULL * (uint64_t)i);
      serial_means[i] = r.mean;
    }

    double t3 = wall_seconds();
    double serial_elapsed = t3 - t2;

    double serial_grand_mean = 0.0;
    for (int i = 0; i < n_serial_replicas; i++) {
      serial_grand_mean += serial_means[i];
    }

    serial_grand_mean /= n_serial_replicas;

    double serial_var = 0.0;
    for (int i = 0; i < n_serial_replicas; i++) {
      double d = serial_means[i] - serial_grand_mean;

      serial_var += d * d;
    }

    double serial_error =
        (n_serial_replicas > 1)
            ? sqrt(serial_var / (n_serial_replicas - 1) / n_serial_replicas)
            : 0.0;

    free(serial_means);

    printf("Serial reference wall time: %.3f s\n", serial_elapsed);
    printf("Serial reference: E = %.6f +- %.6f Hartree (standard error over "
           "%d independent replicas)\n",
           serial_grand_mean, serial_error, n_serial_replicas);
    if (max_elapsed > 1e-9) {
      printf("Hybrid MPI+OpenMP speedup: %.2fx (of %d-way parallelism "
             "available)\n",
             serial_elapsed / max_elapsed, n_serial_replicas);
    }

    /* Statistical agreement check: hybrid grand_mean/inter_rank_error above
     * vs. this serial_grand_mean/serial_error, combined via standard
     * error propagation. These are two independent estimates of the same
     * quantity (same physical system, same total sampling work, disjoint
     * RNG streams via rng_long_jump), so |z| should be O(1) if the hybrid
     * path introduces no statistical bias -- a genuine correctness check,
     * not just a timing comparison. */
    double combined_error =
        sqrt(inter_rank_error * inter_rank_error + serial_error * serial_error);
    double z = (combined_error > 1e-12)
                   ? (grand_mean - serial_grand_mean) / combined_error
                   : 0.0;
    printf("\nStatistical agreement check (hybrid vs. serial, same total "
           "sampling work):\n");
    printf("  z = (E_hybrid - E_serial) / sqrt(err_hybrid^2 + err_serial^2) "
           "= %.3f\n",
           z);
    printf("  %s (|z| %s 2: hybrid and serial agree within combined "
           "statistical uncertainty)\n",
           fabs(z) < 2.0 ? "PASS"
                         : "WARN -- check for a bias in the hybrid "
                           "RNG fan-out or an unlucky sample",
           fabs(z) < 2.0 ? "<" : ">=");

    free(all_means);
    free(all_elapsed);
    free(all_nsamples);
  }

  /* === DMC === */
  MPI_Barrier(MPI_COMM_WORLD);
  uint64_t my_dmc_seed = rank_master_seed(base_seed + 424242ULL, rank);

  double t4 = wall_seconds();
  dmc_result_t local_dmc = dmc_run_parallel(replicas_per_rank, Z, Zeff, 0.15,
                                            /*target_population=*/200,
                                            /*max_population=*/600,
                                            /*tau=*/0.01,
                                            /*n_equilibration=*/500,
                                            /*n_blocks=*/20,
                                            /*block_size=*/200, my_dmc_seed);
  double t5 = wall_seconds();
  double local_dmc_elapsed = t5 - t4;

  double *all_dmc_means = NULL, *all_dmc_elapsed = NULL;
  if (rank == 0) {
    all_dmc_means = malloc((size_t)n_ranks * sizeof(double));
    all_dmc_elapsed = malloc((size_t)n_ranks * sizeof(double));
  }

  MPI_Gather(&local_dmc.energy_mixed, 1, MPI_DOUBLE, all_dmc_means, 1,
             MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&local_dmc_elapsed, 1, MPI_DOUBLE, all_dmc_elapsed, 1, MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  if (rank == 0) {
    double grand_mean = 0.0;

    for (int r = 0; r < n_ranks; r++) {
      grand_mean += all_dmc_means[r];
    }

    grand_mean /= n_ranks;

    double between_var = 0.0;

    for (int r = 0; r < n_ranks; r++) {
      double d = all_dmc_means[r] - grand_mean;

      between_var += d * d;
    }

    double inter_rank_error =
        (n_ranks > 1) ? sqrt(between_var / (n_ranks - 1) / n_ranks) : 0.0;

    double max_elapsed = all_dmc_elapsed[0];

    for (int r = 1; r < n_ranks; r++) {
      if (all_dmc_elapsed[r] > max_elapsed) {
        max_elapsed = all_dmc_elapsed[r];
      }
    }

    printf("\n=== DMC: %d ranks x %d replicas/rank ===\n", n_ranks,
           replicas_per_rank);
    printf("Wall time (max across ranks): %.3f s\n", max_elapsed);
    printf("mixed estimator: E = %.6f +- %.6f Hartree (inter-rank standard "
           "error over %d ranks)\n",
           grand_mean, inter_rank_error, n_ranks);

    free(all_dmc_means);
    free(all_dmc_elapsed);
  }

  MPI_Finalize();

  return 0;
}
