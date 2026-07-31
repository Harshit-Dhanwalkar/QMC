#ifndef QMC_DMC_H
#define QMC_DMC_H

#include "../core/random.h"
#include "vmc.h"
#include <stdint.h>

/*
 * Diffusion Monte Carlo (DMC) for helium ground state, with importance sampling
 * against Slater-Jastrow trial wavefunction (Z'=2 fixed nuclear charge, Jastrow
 * parameter b)
 *
 * All energies in Hartree atomic units
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
} dmc_result_t;

/* Allocate population buffer with given maximum capacity.
 *
 * Returns NULL on allocation failure. count starts at 0; call
 * dmc_population_init to fill it. */
dmc_population_t *dmc_population_alloc(int capacity);
void dmc_population_free(dmc_population_t *pop);

/* Initialize pop with `target_size` walkers (target_size <= pop->capacity
 * required), each drawn via vmc_walker_init. Sets pop->count = target_size. */
void dmc_population_init(dmc_population_t *pop, int target_size, double Zeff,
                         rng_state_t *rng);

/* Drift velocity v = \grad_r(\ln(\Psi_T)) for electron `which` (0 or 1) at
 * walker's current configuration, written into drift_out[3]. */
void dmc_drift_velocity(const vmc_walker_t *w, int which, double Zeff, double b,
                        double drift_out[3]);

/*
 * Attempt single-electron Metropolis-corrected drift-diffusion move:
 * proposes r' = r + \tau*v(r) + \chi, \chi ~ N(0, \tau*I_3), and accepts with
 * probability min(1, |\Psi_T(r')|^2 G(r<-r') / (|\Psi_T(r)|^2 G(r'<-r)))
 *
 * Where
 *   G is drift-diffusion Green's function. Mutates *w on acceptance.
 *
 *Returns 1 if accepted, 0 if rejected.
 */
int dmc_move_electron(vmc_walker_t *w, int which, double Zeff, double b,
                      double tau, rng_state_t *rng);

/*
 * One full DMC generation for single walker: move both electrons once
 * (drift-diffusion Metropolis), then compute branching multiplicity from
 * pre/post-move local energy average and E_T.
 *
 * Returns multiplicity (0 = walker dies, 1 = survives unchanged in count, 2+ =
 * replicated), and writes number of accepted moves (0-2) to *accepted_out if
 * non-NULL. *w is left at its post-move configuration regardless of returned
 * multiplicity.
 */
int dmc_branch_walker(vmc_walker_t *w, double Zeff, double b, double tau,
                      double E_T, rng_state_t *rng, int *accepted_out);

/*
 * Run DMC: initialize population of target_population walkers, equilibrate for
 * n_equilibration generations (discarded), then run n_blocks statistics blocks
 * of block_size generations each, block-averaging both mixed estimator (<E_L>
 * over the population) and growth estimator (E_T history) for reported energy
 * and error.
 *
 * Population is capped at max_population (walkers beyond that are discarded via
 * random subsampling) to bound memory/runtime; max_population should be at
 * least ~3*target_population to avoid frequent, statistics-distorting
 * subsampling under normal branching fluctuations.
 */
dmc_result_t dmc_run(double Zeff, double b, int target_population,
                     int max_population, double tau, int n_equilibration,
                     int n_blocks, int block_size, uint64_t seed);

#endif
