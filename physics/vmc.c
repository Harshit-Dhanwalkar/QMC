#ifndef QMC_VMC_H
#define QMC_VMC_H

#include "../core/random.h"

/*
 * Variational Monte Carlo for the helium ground state.
 * Trial wavefunction (Slater-Jastrow, Z' fixed at bare nuclear charge Z'=2;
 * only Jastrow parameter b is variational):
 *   \Psi_T(r1,r2) = \exp(-Z'(r1+r2)) * \exp( r12 / (2*(1 + b * r12)) )
 * Where
 * r1=|r1_vec|, r2=|r2_vec|, r12=|r1_vec-r2_vec|.
 * The Jastrow prefactor 1/2 is fixed by electron-electron cusp condition; only
 * b is optimized. All energies in Hartree atomic units.
 */

typedef struct {
  double r1[3]; /* electron 1 position */
  double r2[3]; /* electron 2 position */
} vmc_walker_t;

typedef struct {
  double mean;             /* <E_L>, Hartree */
  double error;            /* block-averaged standard error on mean */
  double variance;         /* sample variance of E_L */
  int n_samples;           /* number of retained (post-equilibration) samples */
  double acceptance_rate1; /* electron-1 move acceptance fraction */
  double acceptance_rate2; /* electron-2 move acceptance fraction */
} vmc_result_t;

/* Unnormalized trial wavefunction value at walker configuration. */
double vmc_trial_wavefunction(const vmc_walker_t *w, double Zeff, double b);

/* Local energy E_L = (H \Psi_T)/ \Psi_T at walker configuration. */
double vmc_local_energy(const vmc_walker_t *w, double Zeff, double b);

/* Initialize walker positions by sampling each electron from an exponential
 * radial distribution with rate Z_eff and uniform random direction. */
void vmc_walker_init(vmc_walker_t *w, rng_state_t *rng, double Zeff);

/* Attempt single-electron Metropolis move (which = 0 for electron 1, which = 1
 * for electron 2). Proposes a uniform displacement in [-step_size,
 * step_size]^3, accepts with probability min(1, |\Psi_T(new)|^2 /
 * |\Psi_T(old)|^2). Mutates *w on acceptance.
 * Returns 1 if accepted, 0 if rejected. */
int vmc_metropolis_move_electron(vmc_walker_t *w, int which, double Zeff,
                                 double b, double step_size, rng_state_t *rng);

/* One full sweep = one move attempt per electron. */
void vmc_metropolis_sweep(vmc_walker_t *w, double Zeff, double b,
                          double step_size1, double step_size2,
                          rng_state_t *rng, int *accepted1, int *accepted2);

/* Run VMC: equilibrate for n_equilibration sweeps (discarded), then sample E_L
 * for n_samples sweeps. Error bar computed via block averaging with given
 * block_size to account for serial correlation between successive Metropolis
 * samples.  */
vmc_result_t vmc_run(double Zeff, double b, int n_equilibration, int n_samples,
                     int block_size, double step_size1, double step_size2,
                     uint64_t seed);

/* Optimize b in [b_min, b_max] via golden_section_minimize over
 * vmc_run(...).mean. Uses fixed seed internally per evaluation so objective is
 * deterministic (reproducible minimization), at cost of not re-randomizing
 * across trial b values.
 * Returns optimized mean energy; if b_opt_out is non-NULL, writes the optimal b
 * there. */
double vmc_optimize_b(double Zeff, double b_min, double b_max,
                      int n_equilibration, int n_samples, double step_size1,
                      double step_size2, uint64_t seed, double tol,
                      double *b_opt_out);

#endif
