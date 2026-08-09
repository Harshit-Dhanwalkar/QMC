/*
 * OpenMP-Parallel Quantum Monte Carlo (VMC / DMC / PIMC)
 *
 * vmc_run_parallel / dmc_run_parallel / pimc_run_parallel each run
 * n_replicas fully independent chains -- own walker(s), own equilibration,
 * own statistics -- and combine them. Independence between replicas is
 * exact (not merely statistical): replica i's RNG stream is the master
 * seed's rng_jump()'d i times, a proven non-overlapping 2^128-step advance
 * of the underlying xoshiro256** generator (see core/random.h). The
 * combined error reported is the inter-replica standard error, which does
 * not depend on block_size safely exceeding an unmeasured autocorrelation
 * time the way single-chain block-averaging does.
 *
 * This example runs the same total amount of VMC sampling work serially
 * (n_replicas=1, called n_replicas times in a loop) and in parallel
 * (n_replicas at once via OpenMP) and reports wall-clock time for both, to
 * make the speedup concrete. Build with PLOT_BACKEND=NONE (or any backend)
 * and run with e.g. OMP_NUM_THREADS=4 ./build/eg_39_openmp_qmc to see the
 * effect of thread count -- on a single-core machine (or OMP_NUM_THREADS=1)
 * the two times will be nearly identical, since there is then no actual
 * parallel work being scheduled.
 */

#include "../physics/dmc.h"
#include "../physics/pimc.h"
#include "../physics/vmc.h"
#include <bits/time.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static double wall_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(void) {
  printf(" > OpenMP-Parallel Quantum Monte Carlo\n\n");

#ifdef _OPENMP
  printf("Built with OpenMP. Max threads available: %d\n",
         omp_get_max_threads());
  printf("(Set OMP_NUM_THREADS to control this, e.g. "
         "OMP_NUM_THREADS=4 ./build/eg_39_openmp_qmc)\n\n");
#else
  printf("Built WITHOUT OpenMP (all _parallel functions still run "
         "correctly, just serially).\n\n");
#endif

  /* Helium ground state, same physical parameters throughout. */
  double Z = 2.0, Zeff = 1.6875, b = 0.3;
  int n_equilibration = 1000, n_samples = 60000, block_size = 100;
  double step1 = 0.9, step2 = 0.9;
  int n_replicas = 8;

  printf("=== VMC: %d replicas, %d samples each (%d total) ===\n\n", n_replicas,
         n_samples, n_replicas * n_samples);

  double t0 = wall_seconds();
  double serial_sum_mean = 0.0;
  for (int i = 0; i < n_replicas; i++) {
    /* Each iteration is its own independent single-chain run (n_replicas=1
     * called n_replicas times), same total sampling work as the parallel
     * call below, but executed one chain at a time. */
    vmc_result_t r =
        vmc_run_parallel(1, Z, Zeff, b, n_equilibration, n_samples, block_size,
                         step1, step2, 1000ULL + (uint64_t)i);
    serial_sum_mean += r.mean;
  }
  double t1 = wall_seconds();

  vmc_result_t r_parallel =
      vmc_run_parallel(n_replicas, Z, Zeff, b, n_equilibration, n_samples,
                       block_size, step1, step2, 1000ULL);
  double t2 = wall_seconds();

  printf("Sequential loop of %d single-replica calls: %.3f s "
         "(mean of means = %.6f)\n",
         n_replicas, t1 - t0, serial_sum_mean / n_replicas);
  printf("One vmc_run_parallel(%d, ...) call:          %.3f s "
         "(E = %.6f +- %.6f Hartree)\n",
         n_replicas, t2 - t1, r_parallel.mean, r_parallel.error);
  if (t2 - t1 > 1e-9) {
    printf("Speedup: %.2fx\n\n", (t1 - t0) / (t2 - t1));
  }

  printf("=== DMC: 4 replicas, helium ===\n\n");
  double t3 = wall_seconds();
  dmc_result_t dr = dmc_run_parallel(4, Z, Zeff, 0.15,
                                     /*target_population=*/200,
                                     /*max_population=*/600,
                                     /*tau=*/0.01,
                                     /*n_equilibration=*/500,
                                     /*n_blocks=*/20,
                                     /*block_size=*/200,
                                     /*master_seed=*/2026ULL);
  double t4 = wall_seconds();
  printf("dmc_run_parallel(4, ...): %.3f s, mixed = %.6f +- %.6f Hartree\n\n",
         t4 - t3, dr.energy_mixed, dr.error_mixed);

  printf("=== PIMC: 4 replicas, helium ===\n\n");
  double t5 = wall_seconds();
  int P = 256;
  double tau = 8.0 / P;
  pimc_result_t pr =
      pimc_run_parallel(4, Z, P, tau, /*level=*/4, /*n_equilibration=*/400,
                        /*n_blocks=*/20, /*block_size=*/200,
                        /*master_seed=*/24601ULL);
  double t6 = wall_seconds();
  printf("pimc_run_parallel(4, ...): %.3f s, E = %.6f +- %.6f Hartree\n\n",
         t6 - t5, pr.energy, pr.error);

  printf("Exact helium ground state: -2.9037 Hartree (all three methods "
         "above should bracket or land close to this).\n");

  return 0;
}
