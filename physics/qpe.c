#include "qpe.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "qft.h"
#include "qubits.h"
#include <stdlib.h>

void qpe_run(cvector_t *psi, int n_count, int n_target, const cmatrix_t *U) {
  if (!psi || n_count < 1 || n_target < 1 || !U) {
    return;
  }

  int n_qubits = n_count + n_target;
  long long block = 1LL << n_target;

  if (U->nrows != block || U->ncols != block) {
    return;
  }

  int counting[64];
  for (int k = 0; k < n_count; k++) {
    counting[k] = k;
  }

  int targets[64];
  for (int t = 0; t < n_target; t++) {
    targets[t] = n_count + t;
  }

  for (int k = 0; k < n_count; k++) {
    qstate_apply_gate1(psi, n_qubits, counting[k], hadamard_gate);
  }

  /* Precompute U^(2^i) for i = 0..n_count-1 via repeated squaring:
   * pow[0] = U, pow[i] = pow[i-1] * pow[i - 1]. */
  cmatrix_t **pow = malloc((size_t)n_count * sizeof(cmatrix_t *));
  pow[0] = cmatrix_copy(U);
  for (int i = 1; i < n_count; i++) {
    pow[i] = cmatrix_multiply(pow[i - 1], pow[i - 1]);
  }

  /* Counting qubit k (0 = most - significant / first) controls U^(2^(n_count -
   * 1 - k)), so the first qubit gets the highest power : convention so
   * inverse-QFT readout comes out as a big-endian binary fraction of \phi. */
  for (int k = 0; k < n_count; k++) {
    int exponent = n_count - 1 - k;

    qstate_apply_controlled_unitary(psi, n_qubits, counting[k], targets,
                                    n_target, pow[exponent]);
  }

  for (int i = 0; i < n_count; i++) {
    cmatrix_free(pow[i]);
  }

  free(pow);

  qft_apply_inverse(psi, n_qubits, counting, n_count);
}

qpe_distribution_t *qpe_estimate_phase(int n_count, int n_target,
                                       const cmatrix_t *U,
                                       const complex_t *target_eigenstate) {
  if (n_count < 1 || n_target < 1 || !U || !target_eigenstate) {
    return NULL;
  }

  long long block = 1LL << n_target;
  if (U->nrows != block || U->ncols != block) {
    return NULL;
  }

  int n_qubits = n_count + n_target;
  cvector_t *psi = qstate_alloc(n_qubits);
  if (!psi) {
    return NULL;
  }

  /* Prepare target register in target_eigenstate: overwrite (counting = 0...0,
   * target = *) amplitudes directly, all other amplitudes start at zero from
   * qstate_alloc. */
  for (long long s = 0; s < block; s++) {
    psi->data[s] = target_eigenstate[s];
  }

  qpe_run(psi, n_count, n_target, U);

  qpe_distribution_t *d = malloc(sizeof *d);
  if (!d) {
    cvector_free(psi);

    return NULL;
  }

  long long n_bins = 1LL << n_count;
  d->n_count = n_count;
  d->probabilities = malloc((size_t)n_bins * sizeof(double));
  if (!d->probabilities) {
    free(d);
    cvector_free(psi);

    return NULL;
  }

  /* Marginalize over the target register: \sum |amplitude|^2 over every
   * target-register configuration for each counting-register bitstring. */
  long long dim = 1LL << n_qubits;
  for (long long c = 0; c < n_bins; c++) {
    d->probabilities[c] = 0.0;
  }

  for (long long i = 0; i < dim; i++) {
    long long counting_bits = i >> n_target;

    // NOLINTNEXTLINE(clang-analyzer-core.uninitialized.Assign)
    d->probabilities[counting_bits] += c_abs2(psi->data[i]);
  }

  int best = 0;
  // NOLINTNEXTLINE(clang-analyzer-core.uninitialized.Assign)
  double best_p = d->probabilities[0];
  for (long long c = 1; c < n_bins; c++) {
    if (d->probabilities[c] > best_p) {
      best_p = d->probabilities[c];
      best = (int)c;
    }
  }

  d->phi_estimate = (double)best / (double)n_bins;

  cvector_free(psi);

  return d;
}

void qpe_distribution_free(qpe_distribution_t *d) {
  if (!d) {
    return;
  }

  free(d->probabilities);
  free(d);
}
