#ifndef QMC_QFT_H
#define QMC_QFT_H

#include "../core/vector.h"

/*
 * Quantum Fourier Transform (QFT) and inverse QFT, acting on `n_q` qubits
 * listed in `qubits[]` (MSB-first index convention, same as every qstate_*
 * function in qubits.c) of an n_qubits-qubit state vector `psi`.
 *
 * NOTE: Reference: Quantum Computation and Quantum Information (Nielsen &
 * Chuang sec. 5.1): for each qubit qubits[j] (j=0..n_q-1, most-significant
 * first), apply Hadamard, then a ladder of controlled-phase gates CR_k =
 * diag(1, \exp^{2 * \pi * i / 2^k}) from every later qubit qubits[j+1..n_q-1]
 * acting as control, then a final bit-reversal (swap qubits[j] <->
 * qubits[n_q-1-j]) to match the conventional output bit ordering.
 * qft_apply_inverse runs the same circuit with every phase negated, in reverse
 * gate order (QFT is unitary, so its inverse is its adjoint; for this
 * particular real-diagonal-phase circuit that adjoint is exactly "same gates,
 * phases negated, reversed order").
 *
 * `qubits[]` need not be contiguous or in any particular numeric order : only
 * their given order matters (qubits[0] is treated as the most-significant qubit
 * of the QFT's own local bit convention).
 */
void qft_apply(cvector_t *psi, int n_qubits, const int *qubits, int n_q);
void qft_apply_inverse(cvector_t *psi, int n_qubits, const int *qubits,
                       int n_q);

#endif
