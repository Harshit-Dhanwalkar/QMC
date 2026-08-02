#ifndef QMC_VMC_H
#define QMC_VMC_H

#include "../core/random.h"
#include <stdint.h>

/*
 * Variational Monte Carlo for two-electron atoms/ions (He, H-, Li+, Be2+, ...),
 * Slater-Jastrow ansatz.
 *
 * Nuclear charge Z (in the Hamiltonian's -Z/r1 - Z/r2 potential) and trial
 * orbital exponent Zeff (in wavefunction's \exp(-Zeff * (r1 + r2 )) envelope)
 * are independent parameters: Zeff is variational/screening choice, Z is
 * physical nuclear charge of species being simulated. Only Jastrow parameter b
 * is optimized here; Zeff is fixed by caller.
 *
 * Trial wavefunction:
 *   \Psi_T(r1,r2) = \exp(-Zeff * (r1 + r2)) * \exp(r12 / (2 * (1 + b*r12)))
 *
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
  double variance;         /* sample variance of E_L (raw, not block) */
  int n_samples;           /* number of retained (post-equilibration) samples */
  double acceptance_rate1; /* electron-1 move acceptance fraction */
  double acceptance_rate2; /* electron-2 move acceptance fraction */
} vmc_result_t;

/* Unnormalized trial wavefunction value at walker configuration. */
double vmc_trial_wavefunction(const vmc_walker_t *w, double Zeff, double b);

/* Local energy E_L = (H \Psi_T)/\Psi_T at walker configuration, for
 * two-electron atom/ion of nuclear charge Z with trial orbital exponent Zeff.
 *
 * Returns 0.0 if r1, r2, or r12 is degenerate (< 1e-12), which is
 * probability-zero event during normal sampling
 */
double vmc_local_energy(const vmc_walker_t *w, double Z, double Zeff, double b);

/* Initialize walker positions: each Cartesian component of r1 and r2 drawn from
 * N(0, 1 / Zeff^2).
 */
void vmc_walker_init(vmc_walker_t *w, rng_state_t *rng, double Zeff);

/* Attempt a single-electron Metropolis move (which = 0 for electron 1, which =
 * 1 for electron 2). Proposes uniform displacement in [-step_size, step_size]^3
 * added to that electron's current position, accepts with probability min(1,
 * |\Psi_T(new)|^2 / |\Psi_T(old)|^2). Mutates *w on acceptance.
 *
 * Returns 1 if accepted, 0 if rejected.
 */
int vmc_metropolis_move_electron(vmc_walker_t *w, int which, double Zeff,
                                 double b, double step_size, rng_state_t *rng);

/* One full sweep = one move attempt per electron (electron 1 then electron 2).
 * *accepted1-by-*accepted2 set to 1/0 for this sweep. */
void vmc_metropolis_sweep(vmc_walker_t *w, double Zeff, double b,
                          double step_size1, double step_size2,
                          rng_state_t *rng, int *accepted1, int *accepted2);

/*
 * Run VMC for a two-electron atom/ion of nuclear charge Z, trial orbital
 * exponent Zeff: equilibrate for n_equilibration sweeps (discarded), then
 * sample E_L for n_samples sweeps.
 */
vmc_result_t vmc_run(double Z, double Zeff, double b, int n_equilibration,
                     int n_samples, int block_size, double step_size1,
                     double step_size2, uint64_t seed);

/*
 * Optimize b in [b_min, b_max] via golden_section_minimize over
 * vmc_run(...).mean, for fixed nuclear charge Z and trial orbital exponent
 * Zeff.
 */
double vmc_optimize_b(double Z, double Zeff, double b_min, double b_max,
                      int n_equilibration, int n_samples, double step_size1,
                      double step_size2, uint64_t seed, double tol,
                      double *b_opt_out);

#endif
