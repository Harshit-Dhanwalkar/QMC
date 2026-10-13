/*
 * Multi-target Grover search on a full state vector
 *
 * State vector. Every amplitude is a real double in [0, 1]; complex phases
 * are not needed here because whole algorithm stays in 2D subspace spanned by
 *     |S_marked>   = (1 / \sqrt(M)) \sum_{i in marked} |i>
 *     |S_unmarked> = (1 / \sqrt(N - M)) \sum_{i not in marked} |i>
 * starting from uniform superposition |s> = (|S_marked> + |S_unmarked>) with
 * appropriate coefficients. That subspace is real, and diffusion operator's
 * action on it is a rotation by 2 \theta, so amplitudes never pick up a phase
 * even after many iterations.
 *
 * grover_multi_init:
 *   Writes uniform superposition 1 / \sqrt(N) into amp[0 .. N-1].
 *   Single-precision-avoiding cast `(size_t)1 << n_qubits` keeps state size
 *   computation in size_t, so n_qubits = 12 gives N = 4096 without any
 *   intermediate int overflow. Amplitude itself is computed in double, so
 *   initial norm is 1 to machine precision
 *
 * grover_multi_step (one Grover iteration):
 *   - Phase 1 (oracle) - Every marked amplitude gets its sign flipped in place;
 *                        loop also accumulates sum of amplitudes, which phase 2
 *                        needs
 *   - Phase 2 (diffusion) - Mean of (post-oracle) amplitudes is computed once,
 *                        then each amplitude is reflected about it:
 *                         amp[i] -> 2 * mean - amp[i].
 *   This is 2|s><s| - I expressed in computational basis, and it is whole
 *   inversion-about-mean step. Two passes over array per iteration
 *
 * grover_multi_success:
 *   Sum of |amp[i]|^2 over marked indices. Uses amp[i] * amp[i] rather than
 *   fabs or pow, since amplitudes are real; this is total probability of
 *   measuring a marked item.
 *
 * Analytical references (grover_multi_theta, grover_multi_success_exact,
 * grover_multi_optimal_iterations):
 *   Exact success probability after k iterations is
 *       P_k = \sin^2((2k + 1) \theta)
 *       \sin(theta) = \sqrt(M / N)
 *   and iteration count that comes closest to P = 1 is
 *       k_opt = round(\pi / (4 \theta) - 1/2)
 *   These are evaluated in closed form and are what test suite compares
 *   full-state-vector simulation against
 *
 * Notes on edge cases:
 *   - M = 0 is rejected (\theta would be 0 and algorithm has no fixed point)
 *   - M = N gives \theta = \pi/2, so P_0 = 1 already at k = 0: initial
 *     superposition is entirely on marked items and no iteration is needed.
 *     Closed forms return right value, but optimal-iteration formula
 *     degenerates (k_opt = 0), which is correct
 *   - amp and marked must be disjoint (they are different arrays of different
 *     types in API; two never alias). No input aliasing is possible through
 *     public interface
 */
#include "grover_multi.h"
#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int valid_n(int n_qubits) {
  return n_qubits >= 1 && n_qubits <= GROVER_MAX_QUBITS;
}

int grover_multi_init(double *amp, int n_qubits) {
  if (amp == NULL || !valid_n(n_qubits)) {
    return -1;
  }

  size_t n_items = (size_t)1 << n_qubits;
  double a = 1.0 / sqrt((double)n_items);

  for (size_t i = 0; i < n_items; i++) {
    amp[i] = a;
  }

  return 0;
}

int grover_multi_step(double *amp, int n_qubits, const unsigned char *marked) {
  if (amp == NULL || marked == NULL || !valid_n(n_qubits)) {
    return -1;
  }

  size_t n_items = (size_t)1 << n_qubits;
  double mean = 0.0;

  for (size_t i = 0; i < n_items; i++) {
    if (marked[i]) {
      amp[i] = -amp[i];
    }
    mean += amp[i];
  }

  mean /= (double)n_items;
  for (size_t i = 0; i < n_items; i++) {
    amp[i] = 2.0 * mean - amp[i];
  }

  return 0;
}

double grover_multi_success(const double *amp, int n_qubits,
                            const unsigned char *marked) {
  if (amp == NULL || marked == NULL || !valid_n(n_qubits)) {
    return NAN;
  }

  size_t n_items = (size_t)1 << n_qubits;
  double p = 0.0;

  for (size_t i = 0; i < n_items; i++) {
    if (marked[i]) {
      p += amp[i] * amp[i];
    }
  }

  return p;
}

double grover_multi_theta(int n_qubits, int n_marked) {
  if (!valid_n(n_qubits) || n_marked < 1 || n_marked > (1 << n_qubits)) {
    return NAN;
  }

  return asin(sqrt((double)n_marked / (double)((size_t)1 << n_qubits)));
}

double grover_multi_success_exact(int n_qubits, int n_marked, int k) {
  double theta = grover_multi_theta(n_qubits, n_marked);
  if (isnan(theta) || k < 0) {
    return NAN;
  }

  double s = sin((2.0 * k + 1.0) * theta);

  return s * s;
}

int grover_multi_optimal_iterations(int n_qubits, int n_marked) {
  double theta = grover_multi_theta(n_qubits, n_marked);
  if (isnan(theta)) {
    return -1;
  }

  double k = round(M_PI / (4.0 * theta) - 0.5);

  return k < 0.0 ? 0 : (int)k;
}
