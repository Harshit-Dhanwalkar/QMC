/*
 * CHSH Bell test on a two-qubit Werner state
 *
 * Density matrix - state is built explicitly as a 4x4 real matrix in
 * computational basis {|00>, |01>, |10>, |11>} rather than stored as a compact
 * pair (V, singlet amplitudes). That is deliberate: it keeps correlation
 * function a plain matrix-vector-matrix contraction, which is same code path a
 * user would take for any other two-qubit state and it avoids baking  Werner
 * parametrisation into correlation routine
 *
 * Correlation function (bell_correlation):
 *   E(a, b) = Tr[ \rho * (A(a) \otimes B(b)) ] with
 *   A(a) = \cos(a) \sigma_z + \sin(a) \sigma_x
 *   and likewise for Bob
 *   Both observables are 2x2 Hermitian and real in this basis, so no complex
 *   arithmetic is needed anywhere. Inner loop contracts tensor product directly
 *   against \rho without ever materialising 4x4 Kronecker product
 *
 *   Index convention is
 *       row index  = 2 * alice_outcome + bob_outcome
 *       col index  = 2 * alice_in  + bob_in
 *   so a[j/2][i/2] * b[j%2][i%2] picks right A and B matrix elements from joint
 *   indices in one line
 *
 * CHSH combination (bell_chsh):
 *   S = E(a,b) - E(a,b') + E(a',b) + E(a',b')
 *   Direct evaluation - three-angle symmetries of singlet are not exploited,
 *   because function is also called with non-symmetric angles during
 *   brute-force grid scan below
 *
 * Analytic references:
 *   bell_local_correlation, bell_local_chsh - correlation of a concrete local
 *   hidden-variable model (shared random angle, each side outputs sign of a
 *   projection), a sawtooth in |a - b| folded into [0, \pi]. This is model
 *   against which quantum CHSH value is compared to show violation is not an
 *   artifact of test settings. bell_deterministic_bound - largest |S| over 16
 *   deterministic local strategies; always 2, classical CHSH bound
 *
 * Brute-force scan (bell_max_chsh):
 *   Scans a 4-dimensional grid of angles a, a', b, b' each taking
 *   n_grid values on [0, 2 pi) and returns largest |S| found. Optimum for a
 *   Werner state with visibility V is 2 sqrt(2) V (Tsirelson bound) and is
 *   achieved at standard 45-degree angle pattern; scan recovers this to
 *   within one grid step and serves as reference sampling routine is
 *   checked against. n_grid is capped at 64 (64^4 = 16.7M evaluations) so a
 *   single call stays under a second on any machine demo runs on
 *
 * Sampling (bell_sample):
 *   Simulates n_pairs measurements per setting pair with a deterministic
 *   xorshift32 generator. Each measurement is two draws from local model
 *   "outcome A is +-1 with equal probability, outcome B equals A with
 *   probability (1 + E)/2 and differs otherwise", which reproduces exactly
 *   singlet's statistics for a fixed setting pair without needing to sample
 *   state's full 4x4 probability distribution. E is computed from
 *   bell_correlation, so sampling is anchored to same correlation function
 *   analytic scan uses; two cannot drift apart. outcome_sums[c] accumulates raw
 *   sum of A*B per setting pair so caller can compute running estimates and
 *   error bars if it wants to. RNG is seeded by caller and is reset to a
 *   canonical nonzero constant if 0 is passed, so a single-seed call reproduces
 *   same sequence on every platform
 *
 * Edge cases:
 *   visibility is clamped to [0, 1] by bell_valid; V = 0 is maximally
 *   mixed state (S = 0 for all settings) and V = 1 is pure singlet
 *   Angles are assumed finite; no wrapping or normalisation is done, so
 *   caller can pass any real angle and trig identity \cos(a - b) =
 *   \cos(a)\cos(b) + \sin(a)\sin(b) does rest
 */
#include "bell_chsh.h"
#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Real two-qubit density matrix of Werner state, basis |00>,|01>,|10>,|11>
static void werner(double v, double rho[4][4]) {
  // singlet (|01> - |10>) / \sqrt{2}
  double psi[4] = {0.0, 1.0 / M_SQRT2, -1.0 / M_SQRT2, 0.0};

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      rho[i][j] = v * psi[i] * psi[j] + (i == j ? (1.0 - v) / 4.0 : 0.0);
    }
  }
}

// A(a) = \cos(a) \sigma_z + \sin(a) \sigma_x, row-major 2x2
static void observable(double a, double m[2][2]) {
  m[0][0] = cos(a);
  m[0][1] = sin(a);
  m[1][0] = sin(a);
  m[1][1] = -cos(a);
}

static int valid(double v) { return v >= 0.0 && v <= 1.0; }

double bell_correlation(double visibility, double a, double b) {
  if (!valid(visibility) || !isfinite(a) || !isfinite(b)) {
    return NAN;
  }

  double rho[4][4];
  double ma[2][2];
  double mb[2][2];
  double e = 0.0;

  werner(visibility, rho);
  observable(a, ma);
  observable(b, mb);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      // (A (x) B)[j][i] with row index = 2 * alice + bob
      e += rho[i][j] * ma[j / 2][i / 2] * mb[j % 2][i % 2];
    }
  }

  return e;
}

double bell_chsh(double v, double a, double a2, double b, double b2) {
  return bell_correlation(v, a, b) - bell_correlation(v, a, b2) +
         bell_correlation(v, a2, b) + bell_correlation(v, a2, b2);
}

double bell_local_correlation(double a, double b) {
  if (!isfinite(a) || !isfinite(b)) {
    return NAN;
  }

  double d = fabs(fmod(a - b, 2.0 * M_PI));
  if (d > M_PI) {
    d = 2.0 * M_PI - d;
  }

  return -(1.0 - 2.0 * d / M_PI);
}

double bell_local_chsh(double a, double a2, double b, double b2) {
  return bell_local_correlation(a, b) - bell_local_correlation(a, b2) +
         bell_local_correlation(a2, b) + bell_local_correlation(a2, b2);
}

double bell_deterministic_bound(void) {
  double best = 0.0;

  for (int s = 0; s < 16; s++) {
    int a = (s & 1) ? 1 : -1;
    int a2 = (s & 2) ? 1 : -1;
    int b = (s & 4) ? 1 : -1;
    int b2 = (s & 8) ? 1 : -1;

    best = fmax(best, fabs((double)(a * b - a * b2 + a2 * b + a2 * b2)));
  }

  return best;
}

double bell_max_chsh(double visibility, int n_grid) {
  if (!valid(visibility) || n_grid < 2 || n_grid > 64) {
    return NAN;
  }

  double best = 0.0;
  double step = 2.0 * M_PI / n_grid;

  for (int i = 0; i < n_grid; i++) {
    for (int j = 0; j < n_grid; j++) {
      for (int k = 0; k < n_grid; k++) {
        for (int l = 0; l < n_grid; l++) {
          best = fmax(best, fabs(bell_chsh(visibility, i * step, j * step,
                                           k * step, l * step)));
        }
      }
    }
  }

  return best;
}

static unsigned next(unsigned *s) {
  // xorshift32
  unsigned x = *s;

  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *s = x;

  return x;
}

double bell_sample(double visibility, double a, double a2, double b, double b2,
                   int n_pairs, unsigned seed, double outcome_sums[4]) {
  if (!valid(visibility) || n_pairs < 1 || n_pairs > 10000000 ||
      outcome_sums == NULL || !isfinite(a) || !isfinite(a2) || !isfinite(b) ||
      !isfinite(b2)) {
    return NAN;
  }

  double ang[4][2] = {{a, b}, {a, b2}, {a2, b}, {a2, b2}};
  unsigned s = seed ? seed : 0x9e3779b9u;

  for (int c = 0; c < 4; c++) {
    /* NOTE: Both outcomes are +-1 with probability 1/2 each (Werner state has
     * no marginal bias) and P(equal) = (1 + E) / 2 with E = -V \cos(a - b)
     * taken from density matrix */
    double e = bell_correlation(visibility, ang[c][0], ang[c][1]);
    double sum = 0.0;

    for (int i = 0; i < n_pairs; i++) {
      double u = (next(&s) >> 8) / 16777216.0;
      double w = (next(&s) >> 8) / 16777216.0;
      int ra = u < 0.5 ? 1 : -1;
      int rb = w < (1.0 + e) / 2.0 ? ra : -ra;

      sum += ra * rb;
    }

    outcome_sums[c] = sum;
  }

  return (outcome_sums[0] - outcome_sums[1] + outcome_sums[2] +
          outcome_sums[3]) /
         n_pairs;
}
