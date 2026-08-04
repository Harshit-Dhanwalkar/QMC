/*
Diffusion Monte Carlo for two-electron atoms/ions (He, H-, Li+, Be2+, ...).
*/

#include "dmc.h"
#include "../core/random.h"
#include "vmc.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static void sub3(const double a[3], const double b[3], double out[3]) {
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

static double dot3(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* \ln(\Psi_T)
 * TODO: expose ln_trial_wavefunction from vmc.c/.h instead of duplicating
 * it here, once 2nd consumer (this file) makes duplication */
static double ln_trial_wavefunction(const vmc_walker_t *w, double Zeff,
                                    double b) {
  double r1 = norm3(w->r1);
  double r2 = norm3(w->r2);
  double r12v[3];
  sub3(w->r1, w->r2, r12v);
  double r12 = norm3(r12v);

  return -Zeff * (r1 + r2) + r12 / (2.0 * (1.0 + b * r12));
}

void dmc_drift_velocity(const vmc_walker_t *w, int which, double Zeff, double b,
                        double drift_out[3]) {
  double r1 = norm3(w->r1);
  double r2 = norm3(w->r2);
  double r12v[3];
  sub3(w->r1, w->r2, r12v);
  double s = norm3(r12v);

  if (r1 < 1e-12 || r2 < 1e-12 || s < 1e-12) {
    // configuration guard
    drift_out[0] = drift_out[1] = drift_out[2] = 0.0;

    return;
  }

  double one_plus_bs = 1.0 + b * s;
  double up = 1.0 / (2.0 * one_plus_bs * one_plus_bs); // u'(s)

  if (which == 0) {
    for (int k = 0; k < 3; k++) {
      drift_out[k] = -Zeff * (w->r1[k] / r1) + up * (r12v[k] / s);
    }
  } else {
    for (int k = 0; k < 3; k++) {
      drift_out[k] = -Zeff * (w->r2[k] / r2) - up * (r12v[k] / s);
    }
  }
}

int dmc_move_electron(vmc_walker_t *w, int which, double Zeff, double b,
                      double tau, rng_state_t *rng) {
  if (!w || !rng || (which != 0 && which != 1) || tau <= 0.0) {
    return 0;
  }

  double *moving = (which == 0) ? w->r1 : w->r2;
  double old_pos[3] = {moving[0], moving[1], moving[2]};

  double drift_old[3];
  dmc_drift_velocity(w, which, Zeff, b, drift_old);
  double ln_psi_old = ln_trial_wavefunction(w, Zeff, b);

  double sigma = sqrt(tau);
  double proposed[3];
  for (int k = 0; k < 3; k++) {
    proposed[k] = old_pos[k] + tau * drift_old[k] + sigma * rng_gaussian(rng);
    moving[k] = proposed[k];
  }

  double ln_psi_new = ln_trial_wavefunction(w, Zeff, b);
  double drift_new[3];
  dmc_drift_velocity(w, which, Zeff, b, drift_new);

  /* Green's function ratio for reverse vs forward drift-diffusion move
   * (normalization prefactors (2 * \pi * \tau)^(-3/2) are identical for forward
   * and reverse and cancel in ratio):
   *  fwd = r' - r - \tau * v(r)
   *  bwd = r - r' -
   *  \tau * v(r')
   *  \ln G(r<-r') - \ln G(r'<-r) = (|fwd|^2 - |bwd|^2) / (2 * \tau)
   * Full log acceptance ratio: 2 * (\ln(Psi') - \ln(Psi)) + that
   * Green's-function term. */
  double fwd[3], bwd[3];
  for (int k = 0; k < 3; k++) {
    fwd[k] = proposed[k] - old_pos[k] - tau * drift_old[k];
    bwd[k] = old_pos[k] - proposed[k] - tau * drift_new[k];
  }

  double log_G_ratio = (dot3(fwd, fwd) - dot3(bwd, bwd)) / (2.0 * tau);
  double log_ratio = 2.0 * (ln_psi_new - ln_psi_old) + log_G_ratio;

  int accept;
  if (log_ratio >= 0.0) {
    accept = 1;
  } else {
    double u = rng_uniform(rng);
    accept = (u > 0.0) && (log(u) < log_ratio);
  }

  if (!accept) {
    moving[0] = old_pos[0];
    moving[1] = old_pos[1];
    moving[2] = old_pos[2];
  }

  return accept;
}

int dmc_branch_walker(vmc_walker_t *w, double Z, double Zeff, double b,
                      double tau, double E_T, rng_state_t *rng,
                      int *accepted_out) {
  double E_L_old = vmc_local_energy(w, Z, Zeff, b);

  int a0 = dmc_move_electron(w, 0, Zeff, b, tau, rng);
  int a1 = dmc_move_electron(w, 1, Zeff, b, tau, rng);
  if (accepted_out) {
    *accepted_out = a0 + a1;
  }

  double E_L_new = vmc_local_energy(w, Z, Zeff, b);

  // Trapezoidal average of pre/post-move local energy in branching * weight
  double weight = exp(-tau * (0.5 * (E_L_old + E_L_new) - E_T));

  /* HACK: hard multiplicity cap at 3, following common practical DMC
   * implementations (e.g. Ceperley/Foulkes-style codes).
   * TODO: This is a coarse safeguard, implement smooth energy cutoff (e.g.
   * Umrigar-Nightingale-Runge's velocity/energy cutoff scheme) near
   * singularities instead of flat post-hoc multiplicity clamp. */
  int m = (int)(weight + rng_uniform(rng));
  if (m > 3) {
    m = 3;
  }
  if (m < 0) {
    m = 0;
  }

  return m;
}

dmc_population_t *dmc_population_alloc(int capacity) {
  if (capacity < 1) {
    return NULL;
  }

  dmc_population_t *pop = malloc(sizeof *pop);
  if (!pop) {
    return NULL;
  }

  pop->data = malloc((size_t)capacity * sizeof(vmc_walker_t));
  if (!pop->data) {
    free(pop);
    return NULL;
  }

  pop->count = 0;
  pop->capacity = capacity;

  return pop;
}

void dmc_population_free(dmc_population_t *pop) {
  if (!pop) {
    return;
  }

  free(pop->data);
  free(pop);
}

void dmc_population_init(dmc_population_t *pop, int target_size, double Zeff,
                         rng_state_t *rng) {
  if (!pop || !rng || target_size < 1 || target_size > pop->capacity) {
    return;
  }

  for (int i = 0; i < target_size; i++) {
    vmc_walker_init(&pop->data[i], rng, Zeff);
  }
  pop->count = target_size;
}

/* One full generation: move+branch every walker in `cur`, writing
 survivors/replicas into `next`. Walkers beyond max_population are handled by
 random subsampling back down to target_population.

 Returns population-weighted mean local energy for this generation (for mixed
 estimator and E_T update), and accumulates acceptance counts into
 *accept_sum-by-*move_count.
*/
static double run_one_generation(dmc_population_t *cur, dmc_population_t *next,
                                 double Z, double Zeff, double b, double tau,
                                 double E_T, int target_population,
                                 int max_population, rng_state_t *rng,
                                 long *accept_sum, long *move_count) {
  next->count = 0;
  double E_L_weighted_sum = 0.0;
  long total_copies = 0;

  for (int i = 0; i < cur->count; i++) {
    vmc_walker_t w = cur->data[i];
    int accepted;
    int m = dmc_branch_walker(&w, Z, Zeff, b, tau, E_T, rng, &accepted);

    *accept_sum += accepted;
    *move_count += 2;

    double E_L_new = vmc_local_energy(&w, Z, Zeff, b);
    E_L_weighted_sum += E_L_new * m;
    total_copies += m;

    for (int copy = 0; copy < m; copy++) {
      if (next->count < next->capacity) {
        next->data[next->count] = w;
        next->count++;
      }
      /* NOTE: If next->capacity is exhausted, further copies of walker are
       * silently dropped rather than growing unbounded */
    }
  }

  /*
   * Population control via comb (systematic) resampling back to
   * target_population whenever post-branching count exceeds max_population.
   */
  if (next->count > max_population) {
    int n_pool = next->count;
    double step = (double)n_pool / target_population;
    double offset = rng_uniform(rng) * step;

    vmc_walker_t *resampled =
        malloc((size_t)target_population * sizeof *resampled);
    if (resampled) {
      for (int k = 0; k < target_population; k++) {
        int idx = (int)(offset + k * step);
        if (idx >= n_pool) {
          idx = n_pool - 1; /* floating-point edge-case safety clamp */
        }

        resampled[k] = next->data[idx];
      }

      for (int k = 0; k < target_population; k++) {
        next->data[k] = resampled[k];
      }

      free(resampled);
    }
    /* NOTE: If the allocation failed, next->count is still forced down below so
     * population stays bounded; (rare, transient) walkers left in
     * data[0..target_population-1] from before this block are used as-is rather
     * than leaving population control silently skipped. */

    next->count = target_population;
  }

  return (total_copies > 0) ? E_L_weighted_sum / total_copies : E_T;
}

dmc_result_t dmc_run(double Z, double Zeff, double b, int target_population,
                     int max_population, double tau, int n_equilibration,
                     int n_blocks, int block_size, uint64_t seed) {
  dmc_result_t result = {0};

  if (target_population < 1 || max_population < target_population ||
      tau <= 0.0 || n_blocks < 1 || block_size < 1) {
    return result;
  }

  rng_state_t rng;
  rng_seed(&rng, seed);

  dmc_population_t *pop_a = dmc_population_alloc(max_population);
  dmc_population_t *pop_b = dmc_population_alloc(max_population);
  if (!pop_a || !pop_b) {
    dmc_population_free(pop_a);
    dmc_population_free(pop_b);

    return result;
  }

  dmc_population_init(pop_a, target_population, Zeff, &rng);

  /* E_T feedback gain. kappa~0.1-1 is standard (Umrigar, Nightingale & Runge
   * 1993); 0.1 is chosen here.
   * TODO: expose \kappa as a dmc_run parameter if different (\tau,
   * target_population) regime ever needs different feedback gain to stay stable
   * which is kept internal for now to avoid over-parameterizing single
   * well-tested default. */
  const double kappa = 0.1;
  double E_T = 0.0;
  for (int i = 0; i < pop_a->count; i++) {
    E_T += vmc_local_energy(&pop_a->data[i], Z, Zeff, b);
  }
  E_T /= pop_a->count;

  dmc_population_t *cur = pop_a;
  dmc_population_t *next = pop_b;
  long accept_sum = 0, move_count = 0;

  for (int gen = 0; gen < n_equilibration; gen++) {
    double mean_E_L =
        run_one_generation(cur, next, Z, Zeff, b, tau, E_T, target_population,
                           max_population, &rng, &accept_sum, &move_count);
    int N = next->count;
    E_T = mean_E_L - (kappa / tau) * log((double)N / target_population);

    dmc_population_t *tmp = cur;
    cur = next;
    next = tmp;
  }

  double *block_means_mixed = malloc((size_t)n_blocks * sizeof(double));
  double *block_means_growth = malloc((size_t)n_blocks * sizeof(double));
  if (!block_means_mixed || !block_means_growth) {
    free(block_means_mixed);
    free(block_means_growth);
    dmc_population_free(pop_a);
    dmc_population_free(pop_b);

    return result;
  }

  double pop_size_sum = 0.0;
  int pop_size_count = 0;

  for (int blk = 0; blk < n_blocks; blk++) {
    double sum_mixed = 0.0, sum_growth = 0.0;

    for (int s = 0; s < block_size; s++) {
      double mean_E_L =
          run_one_generation(cur, next, Z, Zeff, b, tau, E_T, target_population,
                             max_population, &rng, &accept_sum, &move_count);
      int N = next->count;
      E_T = mean_E_L - (kappa / tau) * log((double)N / target_population);

      sum_mixed += mean_E_L;
      sum_growth += E_T;
      pop_size_sum += N;
      pop_size_count++;

      dmc_population_t *tmp = cur;
      cur = next;
      next = tmp;
    }

    block_means_mixed[blk] = sum_mixed / block_size;
    block_means_growth[blk] = sum_growth / block_size;
  }

  double mean_mixed = 0.0, mean_growth = 0.0;
  for (int blk = 0; blk < n_blocks; blk++) {
    mean_mixed += block_means_mixed[blk];
    mean_growth += block_means_growth[blk];
  }
  mean_mixed /= n_blocks;
  mean_growth /= n_blocks;

  double err_mixed = 0.0, err_growth = 0.0;
  if (n_blocks > 1) {
    double var_mixed = 0.0, var_growth = 0.0;

    for (int blk = 0; blk < n_blocks; blk++) {
      double dm = block_means_mixed[blk] - mean_mixed;
      double dg = block_means_growth[blk] - mean_growth;
      var_mixed += dm * dm;
      var_growth += dg * dg;
    }

    var_mixed /= (n_blocks - 1);
    var_growth /= (n_blocks - 1);
    err_mixed = sqrt(var_mixed / n_blocks);
    err_growth = sqrt(var_growth / n_blocks);
  }

  free(block_means_mixed);
  free(block_means_growth);

  result.energy_mixed = mean_mixed;
  result.error_mixed = err_mixed;
  result.energy_growth = mean_growth;
  result.error_growth = err_growth;
  result.n_blocks = n_blocks;
  result.mean_population =
      (pop_size_count > 0) ? pop_size_sum / pop_size_count : 0.0;
  result.acceptance_rate =
      (move_count > 0) ? (double)accept_sum / move_count : 0.0;

  dmc_population_free(pop_a);
  dmc_population_free(pop_b);

  return result;
}
