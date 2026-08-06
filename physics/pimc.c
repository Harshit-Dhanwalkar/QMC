/*
Path Integral Monte Carlo for a two-electron atom/ion, Kelbg-regularized Coulomb
pair actions, bisection sampling.
*/

#include "pimc.h"
#include "../core/random.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static int is_power_of_two(int n) { return n > 0 && (n & (n - 1)) == 0; }

double kelbg_potential(double r, double lambda, double q) {
  double x = r / lambda;
  double term1 = (q / r) * (1.0 - exp(-x * x));
  double term2 = q * (sqrt(M_PI) / lambda) * erfc(x);

  return term1 + term2;
}

double kelbg_energy_correction(double r, double lambda, double q) {
  double x = r / lambda;
  double term1 = (q / r) * (1.0 - exp(-x * x));
  double term2 = q * (sqrt(M_PI) / (2.0 * lambda)) * erfc(x);

  return term1 + term2;
}

pimc_walker_t *pimc_walker_alloc(int P) {
  if (P < 1) {
    return NULL;
  }

  pimc_walker_t *w = malloc(sizeof *w);
  if (!w) {
    return NULL;
  }

  w->r1 = malloc((size_t)P * sizeof *w->r1);
  w->r2 = malloc((size_t)P * sizeof *w->r2);
  if (!w->r1 || !w->r2) {
    free(w->r1);
    free(w->r2);
    free(w);

    return NULL;
  }

  w->P = P;

  return w;
}

void pimc_walker_free(pimc_walker_t *w) {
  if (!w) {
    return;
  }

  free(w->r1);
  free(w->r2);
  free(w);
}

void pimc_walker_init(pimc_walker_t *w, rng_state_t *rng, double Z) {
  if (!w || !rng) {
    return;
  }

  double scale = (Z > 0.0) ? 1.0 / Z : 1.0;
  double c1[3] = {0.5, 0.0, 0.0};
  double c2[3] = {-0.5, 0.0, 0.0};

  for (int i = 0; i < w->P; i++) {
    for (int k = 0; k < 3; k++) {
      w->r1[i][k] = c1[k] + rng_gaussian(rng) * scale;
      w->r2[i][k] = c2[k] + rng_gaussian(rng) * scale;
    }
  }
}

/* Potential action over specific set of bead offsets (from `anchor`, modulo P),
 * for given pair of ring paths. Used to compute old-vs-new potential-action
 * difference for bisection segment without touching rest of the ring. `offsets`
 * has `n_offsets` entries.
 */
static double segment_potential_action(double (*r1)[3], double (*r2)[3],
                                       const int *offsets, int n_offsets,
                                       int anchor, int P, double Z,
                                       double tau) {
  double lambda_eN = sqrt(tau);
  double lambda_ee = sqrt(2.0 * tau);

  double S = 0.0;
  for (int k = 0; k < n_offsets; k++) {
    int i = (anchor + offsets[k]) % P;

    double r1n = norm3(r1[i]);
    double r2n = norm3(r2[i]);
    double r12v[3] = {r1[i][0] - r2[i][0], r1[i][1] - r2[i][1],
                      r1[i][2] - r2[i][2]};
    double s = norm3(r12v);

    double v_eN1 = kelbg_potential(r1n, lambda_eN, -Z);
    double v_eN2 = kelbg_potential(r2n, lambda_eN, -Z);
    double v_ee = kelbg_potential(s, lambda_ee, 1.0);

    S += tau * (v_eN1 + v_eN2 + v_ee);
  }

  return S;
}

/* Free-particle bridge midpoint sample: given endpoints x_a, x_b separated by 2
 * * sub_dt in imaginary time (sub_dt = time from either endpoint to midpoint),
 * draws midpoint from N((x_a + x_b)/2, sub_dt/(2 * mass)).
 */
static void bridge_sample(const double xa[3], const double xb[3], double sub_dt,
                          double mass, rng_state_t *rng, double out[3]) {
  double sigma = sqrt(sub_dt / (2.0 * mass));
  for (int k = 0; k < 3; k++) {
    double mean = 0.5 * (xa[k] + xb[k]);

    out[k] = rng_gaussian_scaled(rng, mean, sigma);
  }
}

int pimc_bisection_move(pimc_walker_t *w, int which, double Z, double tau,
                        int level, rng_state_t *rng) {
  if (!w || !rng || (which != 0 && which != 1) || !is_power_of_two(w->P)) {
    return 0;
  }

  int P = w->P;
  int L = 0;
  while ((1 << L) < P) {
    L++;
  } // P = 2^L

  if (level < 1 || level > L) {
    return 0;
  }

  // segment length: seg+1 beads span segment, 2 fixed endpoints + (seg-1)
  // resampled interior
  int seg = 1 << level;

  double (*moving)[3] = (which == 0) ? w->r1 : w->r2;
  double (*other)[3] = (which == 0) ? w->r2 : w->r1;

  int anchor = (int)(rng_uniform(rng) * P);
  if (anchor >= P) {
    anchor = P - 1;
  }

  /* new_seg[0..seg] holds segment's beads at offsets 0..seg from anchor; offset
   * 0 and offset seg are fixed, offsets 1..seg-1 get resampled. */
  double (*new_seg)[3] = malloc((size_t)(seg + 1) * sizeof *new_seg);
  double (*old_seg)[3] = malloc((size_t)(seg + 1) * sizeof *old_seg);
  int *offsets = malloc((size_t)(seg + 1) * sizeof *offsets);
  if (!new_seg || !old_seg || !offsets) {
    free(new_seg);
    free(old_seg);
    free(offsets);

    return 0;
  }

  for (int o = 0; o <= seg; o++) {
    int idx = (anchor + o) % P;
    new_seg[o][0] = old_seg[o][0] = moving[idx][0];
    new_seg[o][1] = old_seg[o][1] = moving[idx][1];
    new_seg[o][2] = old_seg[o][2] = moving[idx][2];
    offsets[o] = o;
  }

  double mass = 1.0; // electron mass, atomic units

  for (int lvl = level; lvl >= 1; lvl--) {
    int step = 1 << lvl;
    int half = step / 2;
    double sub_dt = (double)half * tau;

    for (int start = 0; start < seg; start += step) {
      int end = start + step; // both within [0, seg], no wraparound here
      int mid = start + half;

      bridge_sample(new_seg[start], new_seg[end], sub_dt, mass, rng,
                    new_seg[mid]);
    }
  }

  /* Potential action only needs the OTHER electron's beads at same * offsets
   * (its beads don't move, but the e-e Kelbg term depends on both electrons'
   * positions at each shared bead index). Build a matching "other" segment view
   * (unchanged, same values for old and new).
   */
  double (*other_seg)[3] = malloc((size_t)(seg + 1) * sizeof *other_seg);
  if (!other_seg) {
    free(new_seg);
    free(old_seg);
    free(offsets);

    return 0;
  }

  for (int o = 0; o <= seg; o++) {
    int idx = (anchor + o) % P;
    other_seg[o][0] = other[idx][0];
    other_seg[o][1] = other[idx][1];
    other_seg[o][2] = other[idx][2];
  }

  // Interior offsets only
  int n_interior = seg - 1;
  int *interior_offsets = offsets + 1; // offsets[1..seg-1] == 1..seg-1

  double S_old, S_new;
  if (which == 0) {
    S_old = segment_potential_action(old_seg, other_seg, interior_offsets,
                                     n_interior, 0, seg + 1, Z, tau);
    S_new = segment_potential_action(new_seg, other_seg, interior_offsets,
                                     n_interior, 0, seg + 1, Z, tau);
  } else {
    S_old = segment_potential_action(other_seg, old_seg, interior_offsets,
                                     n_interior, 0, seg + 1, Z, tau);
    S_new = segment_potential_action(other_seg, new_seg, interior_offsets,
                                     n_interior, 0, seg + 1, Z, tau);
  }

  double dS = S_new - S_old;

  int accept;
  if (dS <= 0.0) {
    accept = 1;
  } else {
    double u = rng_uniform(rng);
    accept = (u > 0.0) && (log(u) < -dS);
  }

  if (accept) {
    for (int o = 1; o < seg; o++) {
      int idx = (anchor + o) % P;
      moving[idx][0] = new_seg[o][0];
      moving[idx][1] = new_seg[o][1];
      moving[idx][2] = new_seg[o][2];
    }
  }

  free(new_seg);
  free(old_seg);
  free(offsets);
  free(other_seg);

  return accept;
}

double pimc_energy_estimator(const pimc_walker_t *w, double Z, double tau) {
  if (!w) {
    return 0.0;
  }

  int P = w->P;
  double lambda_eN = sqrt(tau);
  double lambda_ee = sqrt(2.0 * tau);

  const int d = 3;
  const int N = 2; // two electrons

  double kinetic_term = 0.0;
  for (int i = 0; i < P; i++) {
    int ip1 = (i + 1) % P;

    double d1[3] = {w->r1[ip1][0] - w->r1[i][0], w->r1[ip1][1] - w->r1[i][1],
                    w->r1[ip1][2] - w->r1[i][2]};
    double d2[3] = {w->r2[ip1][0] - w->r2[i][0], w->r2[ip1][1] - w->r2[i][1],
                    w->r2[ip1][2] - w->r2[i][2]};

    double sq1 = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2];
    double sq2 = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2];

    kinetic_term += (sq1 + sq2) / (2.0 * tau * tau);
  }

  kinetic_term /= P;

  double potential_term = 0.0;
  for (int i = 0; i < P; i++) {
    double r1n = norm3(w->r1[i]);
    double r2n = norm3(w->r2[i]);
    double r12v[3] = {w->r1[i][0] - w->r2[i][0], w->r1[i][1] - w->r2[i][1],
                      w->r1[i][2] - w->r2[i][2]};
    double s = norm3(r12v);

    potential_term += kelbg_energy_correction(r1n, lambda_eN, -Z);
    potential_term += kelbg_energy_correction(r2n, lambda_eN, -Z);
    potential_term += kelbg_energy_correction(s, lambda_ee, 1.0);
  }

  potential_term /= P;

  /* Leading constant is d * N * P / (2 * \beta),
   *  Which simplifies to d * N / (2 * \tau) since \beta = P * \tau
   */
  return (double)(d * N) / (2.0 * tau) - kinetic_term + potential_term;
}

/*
 * dV_Kelbg/dr at fixed \lambda (\tau held fixed): derived term-by-term from
 *   V_Kelbg(r) = (q/r) * (1 - \exp(-r^2 / \lambda^2))
 *                  + q * (\sqrt(\pi) / \lambda) * erfc(r / \lambda)
 */
static double kelbg_dV_dr(double r, double lambda, double q) {
  double x = r / lambda;

  return -q * (1.0 - exp(-x * x)) / (r * r);
}

double pimc_virial_estimator(const pimc_walker_t *w, double Z, double tau) {
  if (!w) {
    return 0.0;
  }

  int P = w->P;
  double beta = P * tau;
  double lambda_eN = sqrt(tau);
  double lambda_ee = sqrt(2.0 * tau);

  const int d = 3;
  const int N = 2;

  // Each particle's centroid (mean position over its P beads)
  double c1[3] = {0.0, 0.0, 0.0}, c2[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < P; i++) {
    for (int k = 0; k < 3; k++) {
      c1[k] += w->r1[i][k];
      c2[k] += w->r2[i][k];
    }
  }

  for (int k = 0; k < 3; k++) {
    c1[k] /= P;
    c2[k] /= P;
  }

  double sum = 0.0;
  for (int i = 0; i < P; i++) {
    double r1n = norm3(w->r1[i]);
    double r2n = norm3(w->r2[i]);
    double r12v[3] = {w->r1[i][0] - w->r2[i][0], w->r1[i][1] - w->r2[i][1],
                      w->r1[i][2] - w->r2[i][2]};
    double s = norm3(r12v);

    // Same (V + \tau * dV / d\tau) potential term the thermodynamic estimator
    // uses
    double vterm = kelbg_energy_correction(r1n, lambda_eN, -Z) +
                   kelbg_energy_correction(r2n, lambda_eN, -Z) +
                   kelbg_energy_correction(s, lambda_ee, 1.0);

    // grad_1 V and grad_2 V of the total potential at this bead (plain spatial
    // gradient, (\tau / \lambda) held fixed).
    double dV_dr1 = kelbg_dV_dr(r1n, lambda_eN, -Z);
    double dV_ds = kelbg_dV_dr(s, lambda_ee, 1.0);
    double dV_dr2 = kelbg_dV_dr(r2n, lambda_eN, -Z);

    double grad1[3], grad2[3];
    for (int k = 0; k < 3; k++) {
      double r1hat = w->r1[i][k] / r1n;
      double shat = r12v[k] / s;
      grad1[k] = dV_dr1 * r1hat + dV_ds * shat;
    }

    for (int k = 0; k < 3; k++) {
      double r2hat = w->r2[i][k] / r2n;
      double shat = r12v[k] / s; // d(s)/d(r2) = -shat, sign folded in below
      grad2[k] = dV_dr2 * r2hat - dV_ds * shat;
    }

    double gterm = 0.0;
    for (int k = 0; k < 3; k++) {
      gterm += (w->r1[i][k] - c1[k]) * grad1[k];
      gterm += (w->r2[i][k] - c2[k]) * grad2[k];
    }

    gterm *= 0.5;

    sum += vterm + gterm;
  }

  sum /= P;

  return (double)(d * N) / (2.0 * beta) + sum;
}

pimc_result_t pimc_run(double Z, int P, double tau, int level,
                       int n_equilibration, int n_blocks, int block_size,
                       uint64_t seed) {
  pimc_result_t result = {0};

  int L = 0;
  while ((1 << L) < P) {
    L++;
  }

  if (!is_power_of_two(P) || tau <= 0.0 || n_blocks < 1 || block_size < 1 ||
      level < 1 || level > L) {
    return result;
  }

  int seg = 1 << level;
  int moves_per_sweep_per_electron = P / seg; // covers ring once
  if (moves_per_sweep_per_electron < 1) {
    moves_per_sweep_per_electron = 1;
  }

  rng_state_t rng;
  rng_seed(&rng, seed);

  pimc_walker_t *w = pimc_walker_alloc(P);
  if (!w) {
    return result;
  }
  pimc_walker_init(w, &rng, Z);

  long accept_sum = 0, move_count = 0;

  for (int sweep = 0; sweep < n_equilibration; sweep++) {
    for (int m = 0; m < moves_per_sweep_per_electron; m++) {
      accept_sum += pimc_bisection_move(w, 0, Z, tau, level, &rng);
      accept_sum += pimc_bisection_move(w, 1, Z, tau, level, &rng);
      move_count += 2;
    }
  }

  double *block_means = malloc((size_t)n_blocks * sizeof(double));
  double *block_means_virial = malloc((size_t)n_blocks * sizeof(double));
  if (!block_means || !block_means_virial) {
    free(block_means);
    free(block_means_virial);
    pimc_walker_free(w);

    return result;
  }

  for (int blk = 0; blk < n_blocks; blk++) {
    double sum_E = 0.0;
    double sum_E_virial = 0.0;
 
    for (int s = 0; s < block_size; s++) {
      for (int m = 0; m < moves_per_sweep_per_electron; m++) {
        accept_sum += pimc_bisection_move(w, 0, Z, tau, level, &rng);
        accept_sum += pimc_bisection_move(w, 1, Z, tau, level, &rng);
        move_count += 2;
      }
 
      sum_E += pimc_energy_estimator(w, Z, tau);
      sum_E_virial += pimc_virial_estimator(w, Z, tau);
    }
 
    block_means[blk] = sum_E / block_size;
    block_means_virial[blk] = sum_E_virial / block_size;
  }
 
  double mean = 0.0, mean_virial = 0.0;
  for (int blk = 0; blk < n_blocks; blk++) {
    mean += block_means[blk];
    mean_virial += block_means_virial[blk];
  }
  mean /= n_blocks;
  mean_virial /= n_blocks;
 
  double err = 0.0, err_virial = 0.0;
  if (n_blocks > 1) {
    double var = 0.0, var_virial = 0.0;

    for (int blk = 0; blk < n_blocks; blk++) {
      double d = block_means[blk] - mean;
      var += d * d;
      double dv = block_means_virial[blk] - mean_virial;
      var_virial += dv * dv;
    }

    var /= (n_blocks - 1);
    var_virial /= (n_blocks - 1);
    err = sqrt(var / n_blocks);
    err_virial = sqrt(var_virial / n_blocks);
  }
 
  free(block_means);
  free(block_means_virial);
 
  result.energy = mean;
  result.error = err;
  result.energy_virial = mean_virial;
  result.error_virial = err_virial;
  result.n_blocks = n_blocks;
  result.acceptance_rate =
      (move_count > 0) ? (double)accept_sum / move_count : 0.0;
 
  pimc_walker_free(w);
 
  return result;
}
