#ifndef QMC_QPE_H
#define QMC_QPE_H

#include "../core/matrix.h"
#include "../core/vector.h"

/*
 * Quantum Phase Estimation (QPE), Nielsen & Chuang sec. 5.2.
 *
 * NOTE: Circuit: Hadamard every counting qubit, then for k = 0..n_count-1 apply
 * controlled-U^(2^(n_count - 1 - k)) with counting qubit k as control (qubit 0,
 * the first (most-significant) counting qubit, gets the highest power : QPE
 * convention so final inverse-QFT readout comes out directly as a big-endian
 * binary fraction), then inverse-QFT the counting register. Measuring counting
 * register then samples \phi's binary expansion.
 *
 * `\psi` must already be allocated for n_count+n_target qubits with target
 * block prepared in |u> (counting block conventionally starts at |0...0>, i.e.
 * straight from qstate_alloc).
 * NOTE: This function mutates `\psi` in place (does not measure/collapse it,
 * call qstate_measure_qubit or read off qstate_probability afterward, so
 * repeated or partial measurement is left to the caller rather than forced
 * here).
 */
void qpe_run(cvector_t *psi, int n_count, int n_target, const cmatrix_t *U);

typedef struct {
  int n_count;
  double *probabilities; // length 2^n_count: P(measuring counting register
                         // as each of the 2^n_count possible bitstrings)
  double phi_estimate;   // argmax(probabilities) / 2^n_count: the single
                         // most-likely phase estimate
} qpe_distribution_t;

/*
 * HACK: Wrapper: allocate a fresh (n_count+n_target)-qubit state, prepare the
 * target register in the caller-given eigenstate `target_eigenstate` (length
 * 2^n_target, assumed already normalized; an exact eigenstate of U, not just
 * close to one, or the resulting distribution will have real spread beyond the
 * 1/2^n_count expected from n_count's finite precision alone), run qpe_run, and
 * return full counting-register probability distribution plus its argmax phase
 * estimate. Frees the internal state vector; caller owns the returned
 * qpe_distribution_t (free with qpe_distribution_free).
 *
 * Returns NULL on invalid input (n_count<1, n_target<1, U wrong size,
 * target_eigenstate NULL/wrong length).
 */
qpe_distribution_t *qpe_estimate_phase(int n_count, int n_target,
                                       const cmatrix_t *U,
                                       const complex_t *target_eigenstate);

void qpe_distribution_free(qpe_distribution_t *d);

#endif
