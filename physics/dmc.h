#ifndef QMC_DMC_H
#define QMC_DMC_H

#include "../core/random.h"
#include "vmc.h"
#include <stdint.h>

/*
 * Diffusion Monte Carlo (DMC) for two-electron atoms/ions (He, H-, Li+, Be2+,
 * ...), with importance sampling against a Slater-Jastrow trial wavefunction
 *
 * All energies in Hartree atomic units
 *
 * Parallelism: two independent OpenMP layers. dmc_run_parallel parallelizes
 * across n_replicas fully independent DMC runs (embarrassingly parallel, each
 * replica its own rng_jump-derived stream). Additionally, within a single run,
 * dmc_run/dmc_run_with_rng's walker population loop is itself
 * OpenMP-parallelized across walkers each generation. Unlike VMC/PIMC, which
 * are each a single correlated Markov chain (every step depends on the previous
 * one, so there is no legitimate within-chain parallelism opportunity beyond
 * replica level), DMC's population of walkers evolves independently within one
 * generation. Each population slot always uses its own rng_jump-derived stream
 * regardless of which OpenMP thread processes it in any given generation, so
 * results are bit-identical regardless of thread count.
 */

typedef struct {
  vmc_walker_t *data; /* walker buffer, capacity-sized */
  int count;          /* current number of live walkers */
  int capacity;       /* allocated size of data[] */
} dmc_population_t;

typedef struct {
  double energy_mixed;    /* mixed estimator: <E_L> over all walkers,
                             block-averaged */
  double error_mixed;     /* standard error on energy_mixed */
  double energy_growth;   /* growth estimator: block-averaged E_T history */
  double error_growth;    /* standard error on energy_growth */
  int n_blocks;           /* number of statistics blocks actually completed */
  double mean_population; /* average population size over sampling phase */
  double acceptance_rate; /* overall single-electron move acceptance fraction */
  long n_resamples; /* number of generations (across equilibration + * sampling)
                       where post-branching population exceeded max_population
                       and comb resampling engaged */
} dmc_result_t;

/*
 * Allocate population buffer with given maximum capacity.
 *
 * Returns NULL on allocation failure. count starts at 0; call
 * dmc_population_init to fill it.
 */
typedef struct {
  const vmc_walker_t *walker;
  int which;
  double Zeff;
  double b;
} dmc_drift_params_t;

typedef struct {
  int n_blocks;
  int block_size;
} dmc_block_config_t;

dmc_population_t *dmc_population_alloc(int capacity);
void dmc_population_free(dmc_population_t *pop);

/* Initialize pop with `target_size` walkers (target_size <= pop->capacity
 * required), each drawn via vmc_walker_init. Sets pop->count = target_size.
 */
void dmc_population_init(dmc_population_t *pop, int target_size, double Zeff,
                         rng_state_t *rng);

/* Drift velocity v = \grad_r(\ln(\Psi_T)) for electron `which` (0 or 1) at
 * walker's current configuration, written into drift_out[3].
 */
void dmc_drift_velocity(const dmc_drift_params_t *params, double drift_out[3]);

/*
 * Attempt single-electron Metropolis-corrected drift-diffusion move:
 * proposes :
 *     r' = r + \tau * v(r) + \chi, \chi ~ N(0, \tau*I_3),
 * and accepts with probability :
 *    min(1, |\Psi_T(r')|^2 G(r<-r') / (|\Psi_T(r)|^2 G(r'<-r)))
 *
 * Where
 *   G is drift-diffusion Green's function. Mutates *w on acceptance.
 *
 * Returns 1 if accepted, 0 if rejected.
 */
int dmc_move_electron(vmc_walker_t *walker, int which, double Zeff,
                      double param_b, double tau, rng_state_t *rng);

/*
 * One full DMC generation for single walker (nuclear charge Z, trial orbital
 * exponent Zeff): move both electrons once (drift-diffusion Metropolis, using
 * only Zeff), then compute branching multiplicity from pre/post-move local
 * energy average (against true Hamiltonian with nuclear charge Z) and E_T.
 *
 * Returns multiplicity (0 = walker dies, 1 = survives unchanged in count, 2+ =
 * replicated), and writes number of accepted moves (0-2) to *accepted_out if
 * non-NULL. *w is left at its post-move configuration regardless of returned
 * multiplicity.
 */
int dmc_branch_walker(vmc_walker_t *walker, double Z_charge, double Zeff,
                      double param_b, double tau, double E_T, rng_state_t *rng,
                      int *accepted_out);

/*
 * Run DMC for a two-electron atom/ion of nuclear charge Z, trial orbital
 * exponent Zeff: initialize population of target_population walkers,
 * equilibrate for n_equilibration generations (discarded), then run n_blocks
 * statistics blocks of block_size generations each, block-averaging both mixed
 * estimator (<E_L> over population) and growth estimator (E_T history) for
 * reported energy and error.
 *
 * NOTE: Population is capped at max_population (walkers beyond that are
 * discarded via random subsampling) to bound memory/runtime; max_population
 * should be at least ~3*target_population to avoid frequent,
 * statistics-distorting subsampling under normal branching fluctuations.
 */
dmc_result_t dmc_run(double Z_charge, double Zeff, double param_b,
                     int target_population, int max_population, double tau,
                     int n_equilibration, int n_blocks, int block_size,
                     uint64_t seed);

/*
 * Same physics as dmc_run, but runs n_replicas fully independent DMC
 * populations (each with its own walker population, own equilibration, own
 * branching/E_T feedback) and combines them, parallelized over OpenMP threads
 * when built with -fopenmp.
 *
 * NOTE: Replica i's stream is master_seed's rng_jump()'d i times, giving exact
 * independence between replicas. Both error_mixed and error_growth are
 * inter-replica standard error (std of the n_replicas independent replica
 * energies, divided by \sqrt(n_replicas)), which does not depend on
 * n_blocks/block_size being large enough to average out a population's
 * branching-correlation time, unlike dmc_run's single-population block error.
 *
 * Result result.n_blocks = n_blocks * n_replicas (informational total). Falls
 * back to an all-zero result if n_replicas < 1 or the usual dmc_run validity
 * conditions fail, or on allocation failure.
 */
dmc_result_t dmc_run_parallel(int n_replicas, double Z_charge, double Zeff,
                              double param_b, int target_population,
                              int max_population, double tau,
                              int n_equilibration, int n_blocks, int block_size,
                              uint64_t master_seed);

#endif
