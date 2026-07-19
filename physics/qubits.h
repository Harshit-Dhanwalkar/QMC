#ifndef QMC_QUBITS_H
#define QMC_QUBITS_H

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"

/*
 * NOTE: Minimal multi-qubit state vector substrate: single-qubit gates + one
 * entangling two-qubit gate (CNOT) for universal quantum
 * computation (Solovay-Kitaev: arbitrary single-qubit rotations + one
 * entangler generate any N-qubit unitary to arbitrary precision).
 */

// Allocate an n-qubit state initialized to |00...0>, n >= 1 and small enough
// that 2^n doubles fit in memory.
cvector_t *qstate_alloc(int n_qubits);

// P(measuring computational basis state `index`) = |amplitude|^2.
double qstate_probability(const cvector_t *psi, int index);

// Apply an arbitrary single-qubit gate to qubit `target` (0-indexed,
// 0=leftmost) of an n_qubits-qubit state \psi (length 2^n_qubits).
void qstate_apply_gate1(cvector_t *psi, int n_qubits, int target,
                        const complex_t gate[4]);

// Apply gate to `target` only when qubit `control` is |1>: controlled-U
// primitive. control != target required.

// NOTE: CNOT is special case U = \sigma_x (from angular.c).
void qstate_apply_controlled_u(cvector_t *psi, int n_qubits, int control,
                               int target, const complex_t U[4]);

// CNOT(control, target) = qstate_apply_controlled_u with U = \sigma_x.
void qstate_apply_cnot(cvector_t *psi, int n_qubits, int control, int target);

// Hadamard gate: (1/\sqrt2) [[1,1],[1,-1]].
extern const complex_t hadamard_gate[4];

/*
 * Reduced density matrix of single qubit, tracing out all others:
 *   \rho_{ab} = \sum_{rest} conj(\psi[rest,qubit=a]) * \psi[rest,qubit=b]
 * Returns a 2x2 cmatrix_t.
 */
cmatrix_t *qstate_reduced_density_single(const cvector_t *psi, int n_qubits,
                                         int qubit);

/*
 * Von Neumann entropy S = -\sum_i \lambda_i \log2(\lambda_i) of a 2x2
 * Hermitian density matrix, computed via closed-form eigenvalues of 2x2
 * Hermitian matrix (\lambda = (tr +- \sqrt(tr^2 - 4*det)) / 2)
 *
 // NOTE: NOT solved via general eigensolver, so entropy has no dependency on
 // eigen_t's internal.
 *
 * Returns entropy in bits (log base 2): 0 for
 * pure/unentangled qubit, 1 for maximally mixed (maximally entangled) qubit.
 */
double von_neumann_entropy_2x2(cmatrix_t *rho);

#endif
