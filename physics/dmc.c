/*
Diffusion Monte Carlo for two-electron atoms/ions (He, H-, Li+, Be2+, ...).
*/

#include "dmc.h"
#include "../core/random.h"
#include "vmc.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static const double CONST_TWO = 2.0;

static double norm3(const double vec[3]) {
  return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

static void sub3(const double vec_a[3], const double vec_b[3], double out[3]) {
  out[0] = vec_a[0] - vec_b[0];
  out[1] = vec_a[1] - vec_b[1];
  out[2] = vec_a[2] - vec_b[2];
}

static double dot3(const double vec_a[3], const double vec_b[3]) {
  return vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1] + vec_a[2] * vec_b[2];
}

/* \ln(\Psi_T)
 * TODO: expose ln_trial_wavefunction from vmc.c/.h instead of duplicating
 * it here, once 2nd consumer (this file) makes duplication */
static double ln_trial_wavefunction(const vmc_walker_t *walker, double Zeff,
                                    double param_b) {
  double radius1 = norm3(walker->r1);
  double radius2 = norm3(walker->r2);
  double r12v[3];
  sub3(walker->r1, walker->r2, r12v);
  double r12 = norm3(r12v);

  return -Zeff * (radius1 + radius2) +
         r12 / (CONST_TWO * (1.0 + param_b * r12));
}

void dmc_drift_velocity(const dmc_drift_params_t *params, double drift_out[3]) {
  if (!params || !params->walker) {
    return;
  }
  const vmc_walker_t *walker = params->walker;
  double radius1 = norm3(walker->r1);
  double radius2 = norm3(walker->r2);
  double r12v[3];
  sub3(walker->r1, walker->r2, r12v);
  double s = norm3(r12v);

  if (radius1 < 1e-12 || radius2 < 1e-12 || s < 1e-12) {
    // configuration guard
    drift_out[0] = drift_out[1] = drift_out[2] = 0.0;

    return;
  }

  double one_plus_bs = 1.0 + params->b * s;
  double up = 1.0 / (CONST_TWO * one_plus_bs * one_plus_bs); // u'(s)

  if (params->which == 0) {
    for (int k = 0; k < 3; k++) {
      drift_out[k] =
          -params->Zeff * (walker->r1[k] / radius1) + up * (r12v[k] / s);
    }
  } else {
    for (int k = 0; k < 3; k++) {
      drift_out[k] =
          -params->Zeff * (walker->r2[k] / radius2) - up * (r12v[k] / s);
    }
  }
}

int dmc_move_electron(vmc_walker_t *walker, int which, double Zeff,
                      double param_b, double tau, rng_state_t *rng) {
  if (!walker || !rng || (which != 0 && which != 1) || tau <= 0.0) {
    return 0;
  }

  double *moving = (which == 0) ? walker->r1 : walker->r2;
  const double old_pos[3] = {moving[0], moving[1], moving[2]};

  dmc_drift_params_t drift_p_old = {walker, which, Zeff, param_b};
  double drift_old[3];
  dmc_drift_velocity(&drift_p_old, drift_old);
  double ln_psi_old = ln_trial_wavefunction(walker, Zeff, param_b);

  double sigma = sqrt(tau);
  double proposed[3];
  for (int k = 0; k < 3; k++) {
    proposed[k] = old_pos[k] + tau * drift_old[k] + sigma * rng_gaussian(rng);
    moving[k] = proposed[k];
  }

  double ln_psi_new = ln_trial_wavefunction(walker, Zeff, param_b);
  dmc_drift_params_t drift_p_new = {walker, which, Zeff, param_b};
  double drift_new[3];
  dmc_drift_velocity(&drift_p_new, drift_new);

  /*
   * NOTE: Green's function ratio for reverse vs forward drift-diffusion move
   * (normalization prefactors (2 * \pi * \tau)^(-3/2) are identical for forward
   * and reverse and cancel in ratio):
   *  fwd = r' - r - \tau * v(r)
   *  bwd = r - r' - \tau * v(r')
   *  \ln G(r<-r') - \ln G(r'<-r) = (|fwd|^2 - |bwd|^2) / (2 * \tau)
   * Full log acceptance ratio: 2 * (\ln(Psi') - \ln(Psi)) + that
   * Green's-function term.
   */
  double fwd[3], bwd[3];
  for (int k = 0; k < 3; k++) {
    fwd[k] = proposed[k] - old_pos[k] - tau * drift_old[k];
    bwd[k] = old_pos[k] - proposed[k] - tau * drift_new[k];
  }

  double log_G_ratio = (dot3(fwd, fwd) - dot3(bwd, bwd)) / (CONST_TWO * tau);
  double log_ratio = CONST_TWO * (ln_psi_new - ln_psi_old) + log_G_ratio;

  int accept;
  if (log_ratio >= 0.0) {
    accept = 1;
  } else {
    double u_val = rng_uniform(rng);
    accept = (u_val > 0.0) && (log(u_val) < log_ratio);
  }

  if (!accept) {
    moving[0] = old_pos[0];
    moving[1] = old_pos[1];
    moving[2] = old_pos[2];
  }

  return accept;
}

int dmc_branch_walker(vmc_walker_t *walker, double Z_charge, double Zeff,
                      double param_b, double tau, double E_T, rng_state_t *rng,
                      int *accepted_out) {
  double E_L_old = vmc_local_energy(walker, Z_charge, Zeff, param_b);

  int a0 = dmc_move_electron(walker, 0, Zeff, param_b, tau, rng);
  int a1 = dmc_move_electron(walker, 1, Zeff, param_b, tau, rng);
  if (accepted_out) {
    *accepted_out = a0 + a1;
  }

  double E_L_new = vmc_local_energy(walker, Z_charge, Zeff, param_b);

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
 * survivors/replicas into `next`. Walkers beyond max_population are handled by
 * random subsampling back down to target_population.
 *
 * NOTE: Parallelized across walker population. Structured as
 * compute-then-scatter to stay race-free under OpenMP: pass 1 (parallel)
 * evolves every walker independently using its own slot-indexed RNG stream from
 * `walker_streams` (length >= cur->count, provided by caller via rng_jump
 * chaining); pass 2 (serial) does prefix-sum-style scatter into `next->data`
 * and accept/energy reductions, avoiding any race on next->count that
 * concurrent scatter would create. Population-control resampling still uses a
 * single `control_rng` (one rng_uniform call per generation, serial by nature
 * already).
 *
 * Returns population-weighted mean local energy for this generation (for mixed
 * estimator and E_T update), and accumulates acceptance counts into
 * *accept_sum / *move_count.
 */
static double run_one_generation(dmc_population_t *cur, dmc_population_t *next,
                                 double Z_charge, double Zeff, double param_b,
                                 double tau, double E_T, int target_population,
                                 int max_population,
                                 rng_state_t *walker_streams,
                                 rng_state_t *control_rng, long *accept_sum,
                                 long *move_count, long *resample_count) {
  next->count = 0;

  int n = cur->count;
  if (n <= 0) {
    // No walkers to evolve (population died out)
    return E_T;
  }

  vmc_walker_t *evolved = malloc((size_t)n * sizeof(vmc_walker_t));
  int *mult = malloc((size_t)n * sizeof(int));
  int *acc = malloc((size_t)n * sizeof(int));
  double *EL = malloc((size_t)n * sizeof(double));

  if (!evolved || !mult || !acc || !EL) {
    free(evolved);
    free(mult);
    free(acc);
    free(EL);

    return E_T; // allocation failure: fallback with E_T unchanged
  }

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; i++) {
    vmc_walker_t w = cur->data[i];
    int accepted;
    int m = dmc_branch_walker(&w, Z_charge, Zeff, param_b, tau, E_T,
                              &walker_streams[i], &accepted);

    evolved[i] = w;
    mult[i] = m;
    acc[i] = accepted;
    EL[i] = vmc_local_energy(&w, Z_charge, Zeff, param_b);
  }

  double E_L_weighted_sum = 0.0;
  long total_copies = 0;

  for (int i = 0; i < n; i++) {
    *accept_sum += acc[i];
    *move_count += 2;

    E_L_weighted_sum += EL[i] * mult[i];
    total_copies += mult[i];

    for (int copy = 0; copy < mult[i]; copy++) {
      if (next->count < next->capacity) {
        next->data[next->count] = evolved[i];
        next->count++;
      }
      /* NOTE: If next->capacity is exhausted, further copies of walker are
       * silently dropped rather than growing unbounded */
    }
  }

  free(evolved);
  free(mult);
  free(acc);
  free(EL);

  /*
   * Population control via comb (systematic) resampling back to
   * target_population whenever post-branching count exceeds max_population.
   */
  if (next->count > max_population) {
    if (resample_count) {
      (*resample_count)++;
    }

    int n_pool = next->count;
    double step = (double)n_pool / target_population;
    double offset = rng_uniform(control_rng) * step;

    vmc_walker_t *resampled =
        malloc((size_t)target_population * sizeof *resampled);
    if (resampled) {
      for (int k = 0; k < target_population; k++) {
        int idx = (int)(offset + k * step);
        if (idx >= n_pool) {
          idx = n_pool - 1; // floating-point edge-case safety clamp
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

static dmc_result_t dmc_run_with_rng(rng_state_t *rng, double Z_charge,
                                     double Zeff, double param_b,
                                     int target_population, int max_population,
                                     double tau, int n_equilibration,
                                     const dmc_block_config_t *blk_cfg) {
  dmc_result_t result = {0};

  if (!blk_cfg || target_population < 1 || max_population < target_population ||
      tau <= 0.0 || blk_cfg->n_blocks < 1 || blk_cfg->block_size < 1) {
    return result;
  }

  int n_blocks = blk_cfg->n_blocks;
  int block_size = blk_cfg->block_size;

  // Population-control safety valve
  dmc_population_t *pop_a = dmc_population_alloc(3 * max_population);
  dmc_population_t *pop_b = dmc_population_alloc(3 * max_population);
  if (!pop_a || !pop_b) {
    dmc_population_free(pop_a);
    dmc_population_free(pop_b);

    return result;
  }

  dmc_population_init(pop_a, target_population, Zeff, rng);

  /*
   * NOTE: Per-slot RNG streams for the walker population : one independent
   * stream per population slot 0..max_population-1, derived once here via
   * rng_jump chaining from the caller's `rng`, which itself continues to be
   * used serially as the population-control RNG ("control_rng" below) and done
   * before any parallel region opens, so no race deriving these.
   */
  rng_state_t *walker_streams =
      malloc((size_t)max_population * sizeof(rng_state_t));
  if (!walker_streams) {
    dmc_population_free(pop_a);
    dmc_population_free(pop_b);

    return result;
  }

  rng_seed(&walker_streams[0], rng_next_u64(rng));
  for (int i = 1; i < max_population; i++) {
    walker_streams[i] = walker_streams[i - 1];

    rng_jump(&walker_streams[i]);
  }

  /*
   * NOTE: E_T feedback gain. \kappa~0.1-1 is a thereotical value (Reference:
   * Umrigar, Nightingale & Runge 1993); 0.1 is chosen here.
   * TODO: expose \kappa as a dmc_run parameter if different (\tau,
   * target_population) regime ever needs different feedback gain to stay stable
   * which is kept internal for now to avoid over-parameterizing single
   * well-tested default. */
  const double kappa = 0.1;
  double E_T = 0.0;
  for (int i = 0; i < pop_a->count; i++) {
    E_T += vmc_local_energy(&pop_a->data[i], Z_charge, Zeff, param_b);
  }

  E_T /= pop_a->count;

  dmc_population_t *cur = pop_a;
  dmc_population_t *next = pop_b;
  long accept_sum = 0;
  long move_count = 0;
  long resample_count = 0;

  for (int gen = 0; gen < n_equilibration; gen++) {
    double mean_E_L =
        run_one_generation(cur, next, Z_charge, Zeff, param_b, tau, E_T,
                           target_population, max_population, walker_streams,
                           rng, &accept_sum, &move_count, &resample_count);
    int grid_n = next->count;

    E_T = mean_E_L - (kappa / tau) * log((double)grid_n / target_population);

    dmc_population_t *tmp = cur;
    cur = next;
    next = tmp;
  }

  double *block_means_mixed = malloc((size_t)n_blocks * sizeof(double));
  double *block_means_growth = malloc((size_t)n_blocks * sizeof(double));
  if (!block_means_mixed || !block_means_growth) {
    free(block_means_mixed);
    free(block_means_growth);
    free(walker_streams);
    dmc_population_free(pop_a);
    dmc_population_free(pop_b);

    return result;
  }

  double pop_size_sum = 0.0;
  int pop_size_count = 0;

  for (int blk = 0; blk < n_blocks; blk++) {
    double sum_mixed = 0.0;
    double sum_growth = 0.0;

    for (int step = 0; step < block_size; step++) {
      double mean_E_L =
          run_one_generation(cur, next, Z_charge, Zeff, param_b, tau, E_T,
                             target_population, max_population, walker_streams,
                             rng, &accept_sum, &move_count, &resample_count);
      int grid_n = next->count;
      E_T = mean_E_L - (kappa / tau) * log((double)grid_n / target_population);

      sum_mixed += mean_E_L;
      sum_growth += E_T;
      pop_size_sum += grid_n;
      pop_size_count++;

      dmc_population_t *tmp = cur;
      cur = next;
      next = tmp;
    }

    block_means_mixed[blk] = sum_mixed / block_size;
    block_means_growth[blk] = sum_growth / block_size;
  }

  double mean_mixed = 0.0;
  double mean_growth = 0.0;
  for (int blk = 0; blk < n_blocks; blk++) {
    mean_mixed += block_means_mixed[blk];
    mean_growth += block_means_growth[blk];
  }

  mean_mixed /= n_blocks;
  mean_growth /= n_blocks;

  double err_mixed = 0.0;
  double err_growth = 0.0;
  if (n_blocks > 1) {
    double var_mixed = 0.0;
    double var_growth = 0.0;

    for (int blk = 0; blk < n_blocks; blk++) {
      double diff_m = block_means_mixed[blk] - mean_mixed;
      double diff_g = block_means_growth[blk] - mean_growth;

      var_mixed += diff_m * diff_m;
      var_growth += diff_g * diff_g;
    }

    var_mixed /= (n_blocks - 1);
    var_growth /= (n_blocks - 1);
    err_mixed = sqrt(var_mixed / n_blocks);
    err_growth = sqrt(var_growth / n_blocks);
  }

  free(block_means_mixed);
  free(block_means_growth);
  free(walker_streams);

  result.energy_mixed = mean_mixed;
  result.error_mixed = err_mixed;
  result.energy_growth = mean_growth;
  result.error_growth = err_growth;
  result.n_blocks = n_blocks;
  result.mean_population =
      (pop_size_count > 0) ? pop_size_sum / pop_size_count : 0.0;
  result.acceptance_rate =
      (move_count > 0) ? (double)accept_sum / move_count : 0.0;
  result.n_resamples = resample_count;

  dmc_population_free(pop_a);
  dmc_population_free(pop_b);

  return result;
}

dmc_result_t dmc_run(double Z_charge, double Zeff, double param_b,
                     int target_population, int max_population, double tau,
                     int n_equilibration, int n_blocks, int block_size,
                     uint64_t seed) {
  rng_state_t rng;
  rng_seed(&rng, seed);

  dmc_block_config_t blk_cfg = {n_blocks, block_size};

  return dmc_run_with_rng(&rng, Z_charge, Zeff, param_b, target_population,
                          max_population, tau, n_equilibration, &blk_cfg);
}

dmc_result_t dmc_run_parallel(int n_replicas, double Z_charge, double Zeff,
                              double param_b, int target_population,
                              int max_population, double tau,
                              int n_equilibration, int n_blocks, int block_size,
                              uint64_t master_seed) {
  dmc_result_t result = {0};

  if (n_replicas < 1 || target_population < 1 ||
      max_population < target_population || tau <= 0.0 || n_blocks < 1 ||
      block_size < 1) {
    return result;
  }

  // One independent DMC population per replica (one population per MPI
  // rank/thread, statistics combined across independent populations at end)
  rng_state_t *streams = malloc((size_t)n_replicas * sizeof(rng_state_t));
  dmc_result_t *replica_results =
      malloc((size_t)n_replicas * sizeof(dmc_result_t));
  if (!streams || !replica_results) {
    free(streams);
    free(replica_results);

    return result;
  }

  rng_seed(&streams[0], master_seed);
  for (int i = 1; i < n_replicas; i++) {
    streams[i] = streams[i - 1];

    rng_jump(&streams[i]);
  }

  dmc_block_config_t blk_cfg = {n_blocks, block_size};

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n_replicas; i++) {
    replica_results[i] = dmc_run_with_rng(&streams[i], Z_charge, Zeff, param_b,
                                          target_population, max_population,
                                          tau, n_equilibration, &blk_cfg);
  }

  double sum_mixed = 0.0;
  double sum_growth = 0.0;
  double sum_pop = 0.0;
  double sum_acc = 0.0;
  long sum_resamples = 0;
  int valid_replicas = 0;

  for (int i = 0; i < n_replicas; i++) {
    if (replica_results[i].n_blocks <= 0) {
      continue;
    }

    sum_mixed += replica_results[i].energy_mixed;
    sum_growth += replica_results[i].energy_growth;
    sum_pop += replica_results[i].mean_population;
    sum_acc += replica_results[i].acceptance_rate;
    sum_resamples += replica_results[i].n_resamples;

    valid_replicas++;
  }

  if (valid_replicas == 0) {
    free(streams);
    free(replica_results);

    return result;
  }

  double grand_mixed = sum_mixed / valid_replicas;
  double grand_growth = sum_growth / valid_replicas;

  double var_mixed = 0.0;
  double var_growth = 0.0;
  for (int i = 0; i < n_replicas; i++) {
    if (replica_results[i].n_blocks <= 0) {
      continue;
    }

    double diff_m = replica_results[i].energy_mixed - grand_mixed;
    double diff_g = replica_results[i].energy_growth - grand_growth;

    var_mixed += diff_m * diff_m;
    var_growth += diff_g * diff_g;
  }

  double err_mixed = 0.0;
  double err_growth = 0.0;
  if (valid_replicas > 1) {
    var_mixed /= (valid_replicas - 1);
    var_growth /= (valid_replicas - 1);
    err_mixed = sqrt(var_mixed / valid_replicas);
    err_growth = sqrt(var_growth / valid_replicas);
  }

  result.energy_mixed = grand_mixed;
  result.error_mixed = err_mixed;
  result.energy_growth = grand_growth;
  result.error_growth = err_growth;
  result.n_blocks = n_blocks * valid_replicas;
  result.mean_population = sum_pop / valid_replicas;
  result.acceptance_rate = sum_acc / valid_replicas;
  result.n_resamples = sum_resamples;

  free(streams);
  free(replica_results);

  return result;
}
