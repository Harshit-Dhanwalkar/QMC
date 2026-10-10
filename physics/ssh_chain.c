/*
 * Su-Schrieffer-Heeger chain
 *
 * Bond table (ssh_chain_bonds):
 *   For an N-cell chain with two sites per cell there are 2N - 1 bonds. Bond i
 *   is v if i is even, w if i is odd (alternating strong/weak), each multiplied
 *   by an independent factor (1 + noise * u) with u uniform in [-1, 1). noise =
 *   0 gives clean dimerized chain used by topological phase diagram; small
 *   noise is a controlled way to break exact chiral symmetry without destroying
 *   invariant
 *
 * Real-space spectrum (ssh_chain_solve):
 *   Builds 2N x 2N tridiagonal Hamiltonian with on-site +m / -m (sublattice
 *   staggering, ssh_chain_onsite) and off-diagonal bonds from ssh_chain_bonds,
 *   then calls tridiag_eigh (or tridiag_eigvals if no eigenvectors are
 *   requested). Eigenvalues are ascending; eigenvectors live column-major
 *   inside returned eigen_t and, when requested, are copied out into caller's
 *   row-major buffer as vectors[s * 2N + i] = Re(<i | s>). Real-space and
 *   analytic-parameter inputs are validated before solver runs so solver itself
 *   sees only finite numbers
 *
 * Analytic references:
 *   ssh_bulk_energy        - \sqrt(v^2 + w^2 + 2 v w \cos(k)), two-band bulk
 *                            dispersion for clean chain
 *   ssh_winding_number     - winding of (v + w \exp(ik) around origin by
 *                            unwrapping \arg(v + w \exp(ik) over 4096 evenly
 *                            spaced k. Returns ±1 in topological phase and 0 in
 *                            trivial phase; -1 for |v| == |w| where winding is
 *                            undefined
 *   ssh_zak_phase          - many-body Zak phase of lower band computed as \arg
 *                            \prod_k <u_k|u_{k + 1}> using explicit lower-band
 *                            spinor (1, -\exp(-i\phi_k) / \sqrt(2). Returns 0
 *                            or \pi (two allowed values by inversion symmetry)
 *                            and NaN at gap closure
 *   ssh_edge_decay_length  - 1 / \log(|w|/|v|), decay length of topological
 *                            edge state in |w| > |v| phase
 */
#include "ssh_chain.h"
#include "../core/linalg/tridiag_eigh.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static unsigned long long next_u64(unsigned long long *state) {
  unsigned long long z = (*state += 0x9E3779B97F4A7C15ULL);

  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;

  return z ^ (z >> 31);
}

int ssh_chain_bonds(double *bond, int n_cells, double v, double w, double noise,
                    unsigned long long seed) {
  if (!bond || n_cells < 1 || !isfinite(v) || !isfinite(w) ||
      !(noise >= 0.0 && noise < 1.0)) {
    return -1;
  }

  unsigned long long state = seed;
  for (int i = 0; i < 2 * n_cells - 1; i++) {
    double u =
        2.0 * ((double)(next_u64(&state) >> 11) * (1.0 / 9007199254740992.0)) -
        1.0;

    bond[i] = ((i % 2 == 0) ? v : w) * (1.0 + noise * u);
  }

  return 0;
}

int ssh_chain_onsite(double *eps, int n_cells, double m) {
  if (!eps || n_cells < 1 || !isfinite(m)) {
    return -1;
  }

  for (int i = 0; i < 2 * n_cells; i++) {
    eps[i] = (i % 2 == 0) ? m : -m;
  }

  return 0;
}

int ssh_chain_solve(int n_cells, const double *bond, const double *eps,
                    double *energies, double *vectors) {
  if (!bond || !energies || n_cells < 1) {
    return -1;
  }

  int n = 2 * n_cells;
  double *diag = calloc((size_t)n, sizeof(double));

  if (!diag) {
    return -2;
  }

  for (int i = 0; i < n - 1; i++) {
    if (!isfinite(bond[i])) {
      free(diag);

      return -1;
    }
  }

  for (int i = 0; i < n; i++) {
    diag[i] = eps ? eps[i] : 0.0;
    if (!isfinite(diag[i])) {
      free(diag);

      return -1;
    }
  }

  eigen_t *eig =
      vectors ? tridiag_eigh(diag, bond, n) : tridiag_eigvals(diag, bond, n);

  free(diag);

  if (!eig) {
    return -2;
  }

  for (int s = 0; s < n; s++) {
    energies[s] = eig->eigenvalues[s];
    if (vectors) {
      for (int i = 0; i < n; i++) {
        vectors[(size_t)s * n + i] = eig->eigenvectors->data[i * n + s].re;
      }
    }
  }

  eigen_free(eig);

  return 0;
}

double ssh_bulk_energy(double v, double w, double k) {
  if (!isfinite(v) || !isfinite(w) || !isfinite(k)) {
    return NAN;
  }

  double e2 = v * v + w * w + 2.0 * v * w * cos(k);

  return sqrt(e2 > 0.0 ? e2 : 0.0);
}

int ssh_winding_number(double v, double w) {
  if (!isfinite(v) || !isfinite(w) || fabs(fabs(v) - fabs(w)) < 1e-12) {
    return -1;
  }

  const int nk = 4096;
  double total = 0.0;
  double prev = atan2(w * 0.0, v + w); // k = 0: h = v + w

  for (int j = 1; j <= nk; j++) {
    double k = 2.0 * M_PI * j / nk;
    double phi = atan2(w * sin(k), v + w * cos(k));
    double d = phi - prev;

    while (d > M_PI) {
      d -= 2.0 * M_PI;
    }
    while (d < -M_PI) {
      d += 2.0 * M_PI;
    }

    total += d;
    prev = phi;
  }

  return (int)lround(total / (2.0 * M_PI));
}

double ssh_zak_phase(double v, double w) {
  if (!isfinite(v) || !isfinite(w) || (v == 0.0 && w == 0.0)) {
    return NAN;
  }

  const int nk = 4096;
  double re = 1.0;
  double im = 0.0;

  /* lower-band eigenvector u_k = (1, -\exp(-i \phi_k) / \sqrt(2))
   * overlap of neighbouring k points = (1 + \exp(i (\phi_k - \phi_k'))) / 2 */
  for (int j = 0; j < nk; j++) {
    double k0 = 2.0 * M_PI * j / nk;
    double k1 = 2.0 * M_PI * (j + 1) / nk;
    double p0 = atan2(w * sin(k0), v + w * cos(k0));
    double p1 = atan2(w * sin(k1), v + w * cos(k1));
    double ore = 0.5 * (1.0 + cos(p0 - p1));
    double oim = 0.5 * sin(p0 - p1);
    double nre = re * ore - im * oim;
    double nim = re * oim + im * ore;

    re = nre;
    im = nim;
  }

  return fabs(atan2(-im, re));
}

double ssh_edge_decay_length(double v, double w) {
  if (!isfinite(v) || !isfinite(w) || !(fabs(v) > 0.0) ||
      !(fabs(w) > fabs(v))) {
    return NAN;
  }

  return 1.0 / log(fabs(w) / fabs(v));
}
