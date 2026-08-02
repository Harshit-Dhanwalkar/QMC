/*
Variational Monte Carlo for helium ground state.
*/

#include "vmc.h"
#include "../core/random.h"
#include "variational.h"
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

// ln(\Psi_T) = -Zeff * (r1 + r2) + r12/(2 * (1 + b * r12))
static double ln_trial_wavefunction(const vmc_walker_t *w, double Zeff,
                                    double b) {
  double r1 = norm3(w->r1);
  double r2 = norm3(w->r2);
  double r12v[3];
  sub3(w->r1, w->r2, r12v);
  double r12 = norm3(r12v);

  return -Zeff * (r1 + r2) + r12 / (2.0 * (1.0 + b * r12));
}

double vmc_trial_wavefunction(const vmc_walker_t *w, double Zeff, double b) {
  if (!w) {
    return 0.0;
  }

  return exp(ln_trial_wavefunction(w, Zeff, b));
}

/*
 * Derivation summary (validated against a finite-difference Laplacian to
 * ~1e-5 for Z=Zeff [regression] and Z != Zeff [H-, Li+, Be2+]):
 *
 * NOTE: Let f = \ln(\Psi_T) = -Zeff * (r1 + r2) + u(s), s = |r1-r2|, u(s) =
 * s/(2 * (1 + bs)). For each electron i, -1/2 * lap_i(\Psi) / \Psi = (-1/2) *
 * (lap_i f + |\grad_i f|^2). Summing i=1,2 and adding potential (-Z/r1 - Z/r2 +
 * 1/s), the Zeff^2 orbital-kinetic term and the (-Z/r1 - Z/r2) potential term
 * combine to leave a residual (Zeff - Z) * (1/r1 + 1/r2) cross term whenever
 * Zeff != Z. Every other term is independent of Z, since only the orbital
 * exponent Zeff enters wavefunction envelope.
 */
double vmc_local_energy(const vmc_walker_t *w, double Z, double Zeff,
                        double b) {
  if (!w) {
    return 0.0;
  }

  double r1 = norm3(w->r1);
  double r2 = norm3(w->r2);
  double r12v[3];
  sub3(w->r1, w->r2, r12v);
  double s = norm3(r12v); // r12

  if (r1 < 1e-12 || r2 < 1e-12 || s < 1e-12) {
    // Probability-Zero
    return 0.0; // degenerate configuration
  }

  double r1hat[3] = {w->r1[0] / r1, w->r1[1] / r1, w->r1[2] / r1};
  double r2hat[3] = {w->r2[0] / r2, w->r2[1] / r2, w->r2[2] / r2};
  double r12hat[3] = {r12v[0] / s, r12v[1] / s, r12v[2] / s};
  double dot_diff[3] = {r1hat[0] - r2hat[0], r1hat[1] - r2hat[1],
                        r1hat[2] - r2hat[2]};
  double angular = dot3(dot_diff, r12hat);

  double one_plus_bs = 1.0 + b * s;

  double E_L = -Zeff * Zeff;
  E_L += (Zeff - Z) * (1.0 / r1 + 1.0 / r2);
  E_L += b / (one_plus_bs * one_plus_bs * one_plus_bs);
  E_L += -1.0 / (s * one_plus_bs * one_plus_bs);
  E_L += -1.0 / (4.0 * one_plus_bs * one_plus_bs * one_plus_bs * one_plus_bs);
  E_L += (Zeff / (2.0 * one_plus_bs * one_plus_bs)) * angular;
  E_L += 1.0 / s;

  return E_L;
}

void vmc_walker_init(vmc_walker_t *w, rng_state_t *rng, double Zeff) {
  if (!w || !rng || Zeff <= 0.0) {
    return;
  }

  double scale = 1.0 / Zeff;
  for (int k = 0; k < 3; k++) {
    w->r1[k] = rng_gaussian(rng) * scale;
    w->r2[k] = rng_gaussian(rng) * scale;
  }
}

int vmc_metropolis_move_electron(vmc_walker_t *w, int which, double Zeff,
                                 double b, double step_size, rng_state_t *rng) {
  if (!w || !rng || (which != 0 && which != 1)) {
    return 0;
  }

  double *moving = (which == 0) ? w->r1 : w->r2;
  double old_pos[3] = {moving[0], moving[1], moving[2]};
  double old_ln = ln_trial_wavefunction(w, Zeff, b);

  for (int k = 0; k < 3; k++) {
    moving[k] += rng_uniform_range(rng, -step_size, step_size);
  }

  double new_ln = ln_trial_wavefunction(w, Zeff, b);
  double log_ratio = 2.0 * (new_ln - old_ln); // |\Psi_new / \Psi_old|^2

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

void vmc_metropolis_sweep(vmc_walker_t *w, double Zeff, double b,
                          double step_size1, double step_size2,
                          rng_state_t *rng, int *accepted1, int *accepted2) {
  int a1 = vmc_metropolis_move_electron(w, 0, Zeff, b, step_size1, rng);
  int a2 = vmc_metropolis_move_electron(w, 1, Zeff, b, step_size2, rng);

  if (accepted1) {
    *accepted1 = a1;
  }
  if (accepted2) {
    *accepted2 = a2;
  }
}

vmc_result_t vmc_run(double Z, double Zeff, double b, int n_equilibration,
                     int n_samples, int block_size, double step_size1,
                     double step_size2, uint64_t seed) {
  vmc_result_t result = {0};

  if (n_samples <= 0 || block_size <= 0) {
    return result;
  }

  rng_state_t rng;
  rng_seed(&rng, seed);

  vmc_walker_t w;
  vmc_walker_init(&w, &rng, Zeff);

  for (int i = 0; i < n_equilibration; i++) {
    int a1, a2;
    vmc_metropolis_sweep(&w, Zeff, b, step_size1, step_size2, &rng, &a1, &a2);
  }

  int n_blocks = n_samples / block_size;
  double *block_means = NULL;
  if (n_blocks > 0) {
    block_means = malloc((size_t)n_blocks * sizeof(double));
    if (!block_means) {
      return result;
    }
  }

  double sum_E = 0.0, sum_E2 = 0.0;
  int sample_count = 0;
  int acc1_total = 0, acc2_total = 0;

  for (int blk = 0; blk < n_blocks; blk++) {
    double block_sum = 0.0;

    for (int s = 0; s < block_size; s++) {
      int a1, a2;
      vmc_metropolis_sweep(&w, Zeff, b, step_size1, step_size2, &rng, &a1, &a2);
      acc1_total += a1;
      acc2_total += a2;

      double E = vmc_local_energy(&w, Z, Zeff, b);
      block_sum += E;
      sum_E += E;
      sum_E2 += E * E;
      sample_count++;
    }
    block_means[blk] = block_sum / block_size;
  }

  double mean = (sample_count > 0) ? sum_E / sample_count : 0.0;
  double variance =
      (sample_count > 0) ? (sum_E2 / sample_count - mean * mean) : 0.0;

  double error = 0.0;
  if (n_blocks > 1) {
    double block_var = 0.0;

    for (int blk = 0; blk < n_blocks; blk++) {
      double d = block_means[blk] - mean;
      block_var += d * d;
    }

    block_var /= (n_blocks - 1);
    error = sqrt(block_var / n_blocks);
  }

  free(block_means);

  result.mean = mean;
  result.error = error;
  result.variance = variance;
  result.n_samples = sample_count;
  result.acceptance_rate1 =
      (sample_count > 0) ? (double)acc1_total / sample_count : 0.0;
  result.acceptance_rate2 =
      (sample_count > 0) ? (double)acc2_total / sample_count : 0.0;

  return result;
}

typedef struct {
  double Z;
  double Zeff;
  int n_equilibration;
  int n_samples;
  int block_size;
  double step_size1;
  double step_size2;
  uint64_t seed;
} vmc_optimize_closure_t;

static double vmc_optimize_eval(double b, void *params) {
  vmc_optimize_closure_t *c = (vmc_optimize_closure_t *)params;
  vmc_result_t r =
      vmc_run(c->Z, c->Zeff, b, c->n_equilibration, c->n_samples, c->block_size,
              c->step_size1, c->step_size2, c->seed);

  return r.mean;
}

double vmc_optimize_b(double Z, double Zeff, double b_min, double b_max,
                      int n_equilibration, int n_samples, double step_size1,
                      double step_size2, uint64_t seed, double tol,
                      double *b_opt_out) {
  const int block_size = 200;

  vmc_optimize_closure_t closure = {Z,          Zeff,       n_equilibration,
                                    n_samples,  block_size, step_size1,
                                    step_size2, seed};

  double b_opt =
      golden_section_minimize(b_min, b_max, vmc_optimize_eval, &closure, tol);

  if (b_opt_out) {
    *b_opt_out = b_opt;
  }

  vmc_result_t final = vmc_run(Z, Zeff, b_opt, n_equilibration, n_samples,
                               block_size, step_size1, step_size2, seed);

  return final.mean;
}
