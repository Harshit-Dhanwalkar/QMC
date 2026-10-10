/*
 * qmc_ssh.c - WebAssembly front end for SSH topological-chain demo
 *
 * Physics: physics/ssh_chain.c (tridiagonal diagonalisation of chain, bulk
 * band, winding number, Zak phase). This file owns one chain of SSH_CELLS unit
 * cells with intra-cell hopping v = 1 and inter-cell hopping w, and buffers
 * JavaScript reads
 */

#include "../physics/ssh_chain.h"
#include <math.h>
#include <stddef.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { SSH_CELLS = 20, SSH_SITES = 2 * SSH_CELLS, SSH_SWEEP_MAX = 160 };

static double g_bond[SSH_SITES - 1];
static double g_eps[SSH_SITES];
static double g_energy[SSH_SITES];
static double g_vec[SSH_SITES * SSH_SITES];
static double g_sweep[SSH_SWEEP_MAX * SSH_SITES];
static double g_w = 1.0;
static double g_noise;
static double g_m;
static int g_seed = 1;

QMC_EXPORT(qmc_ssh_cells)
int qmc_ssh_cells(void) { return SSH_CELLS; }

QMC_EXPORT(qmc_ssh_sweep_max)
int qmc_ssh_sweep_max(void) { return SSH_SWEEP_MAX; }

// Ascending eigenvalues of current chain (2 * cells)
QMC_EXPORT(qmc_ssh_energies)
double *qmc_ssh_energies(void) { return g_energy; }

// Eigenvector s on site i at [s * sites + i]
QMC_EXPORT(qmc_ssh_vectors)
double *qmc_ssh_vectors(void) { return g_vec; }

// 2 * cells - 1 bond strengths (noise included)
QMC_EXPORT(qmc_ssh_bonds)
double *qmc_ssh_bonds(void) { return g_bond; }

/* Energy of sweep point j, state s at [j * sites + s]. */
QMC_EXPORT(qmc_ssh_sweep_buffer)
double *qmc_ssh_sweep_buffer(void) { return g_sweep; }

/* Build and diagonalise chain: v = 1, inter-cell hopping w in [0, 3],
 * hopping disorder `noise` in [0, 1), staggered potential m in [-1, 1] and a
 * seed for disorder
 *
 * Returns 0 on success */
QMC_EXPORT(qmc_ssh_set)
int qmc_ssh_set(double w, double noise, double m, int seed) {
  if (!(w >= 0.0 && w <= 3.0) || !(noise >= 0.0 && noise < 1.0) ||
      !(m >= -1.0 && m <= 1.0)) {
    return -1;
  }

  if (ssh_chain_bonds(g_bond, SSH_CELLS, 1.0, w, noise,
                      (unsigned long long)(unsigned)seed) != 0 ||
      ssh_chain_onsite(g_eps, SSH_CELLS, m) != 0 ||
      ssh_chain_solve(SSH_CELLS, g_bond, g_eps, g_energy, g_vec) != 0) {
    return -1;
  }

  g_w = w;
  g_noise = noise;
  g_m = m;
  g_seed = seed;

  return 0;
}

/* Spectrum of `count` chains with w from w0 to w1, into sweep buffer
 *
 * Returns 0 on success */
QMC_EXPORT(qmc_ssh_sweep)
int qmc_ssh_sweep(double w0, double w1, int count) {
  double bond[SSH_SITES - 1];

  if (count < 2 || count > SSH_SWEEP_MAX || !(w0 >= 0.0) || !(w1 <= 3.0) ||
      !(w1 > w0)) {
    return -1;
  }

  for (int j = 0; j < count; j++) {
    double w = w0 + (w1 - w0) * j / (count - 1);

    if (ssh_chain_bonds(bond, SSH_CELLS, 1.0, w, g_noise,
                        (unsigned long long)(unsigned)g_seed) != 0 ||
        ssh_chain_solve(SSH_CELLS, bond, g_eps, g_sweep + (size_t)j * SSH_SITES,
                        NULL) != 0) {
      return -1;
    }
  }

  return 0;
}

QMC_EXPORT(qmc_ssh_winding)
int qmc_ssh_winding(void) { return ssh_winding_number(1.0, g_w); }

QMC_EXPORT(qmc_ssh_zak)
double qmc_ssh_zak(void) { return ssh_zak_phase(1.0, g_w); }

// Edge-state decay length in cells (NaN in trivial phase)
QMC_EXPORT(qmc_ssh_xi)
double qmc_ssh_xi(void) { return ssh_edge_decay_length(1.0, g_w); }

// Bulk gap 2 |v - w| of infinite chain
QMC_EXPORT(qmc_ssh_gap)
double qmc_ssh_gap(void) { return 2.0 * fabs(1.0 - g_w); }

// Upper bulk band at wave number k for current w
QMC_EXPORT(qmc_ssh_bulk)
double qmc_ssh_bulk(double k) { return ssh_bulk_energy(1.0, g_w, k); }

// Probability of state s within the outermost `edge` cells at either end
QMC_EXPORT(qmc_ssh_edge_weight)
double qmc_ssh_edge_weight(int s, int edge) {
  double sum = 0.0;

  if (s < 0 || s >= SSH_SITES || edge < 1 || 2 * edge > SSH_CELLS) {
    return -1.0;
  }

  for (int i = 0; i < 2 * edge; i++) {
    double a = g_vec[s * SSH_SITES + i];
    double b = g_vec[s * SSH_SITES + SSH_SITES - 1 - i];

    sum += a * a + b * b;
  }

  return sum;
}
