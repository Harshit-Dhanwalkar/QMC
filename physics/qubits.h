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
 *
 * Returns a 2x2 cmatrix_t.
 */
cmatrix_t *qstate_reduced_density_single(const cvector_t *psi, int n_qubits,
                                         int qubit);

/*
 * Von Neumann entropy S = -\sum_i \lambda_i \log2(\lambda_i) of a 2x2
 * Hermitian density matrix, computed via closed-form eigenvalues of 2x2
 * Hermitian matrix (\lambda = (tr \pm \sqrt(tr^2 - 4 * det)) / 2)
 *
 // NOTE: NOT solved via general eigensolver, so entropy has no dependency on
 // eigen_t's internal.
 *
 * Returns entropy in bits (log base 2): 0 for
 * pure/unentangled qubit, 1 for maximally mixed (maximally entangled) qubit.
 */
double von_neumann_entropy_2x2(cmatrix_t *rho);

/*
 * Projective measurement in computational basis, on pure state vector.
 * Samples an outcome index in [0, \psi->n) with probability |\psi[outcome]|^2,
 * using caller-supplied uniform random number u in [0,1) so RNG/seeding stays
 * under caller control. Collapses \psi in place to pure basis state |outcome>
 * (\psi[outcome] = 1, every other amplitude = 0).
 *
 * Returns sampled outcome index, or -1 on invalid input.
 */
int qstate_measure(cvector_t *psi, double u);

/*
 * Projective measurement of single qubit `target` (0-indexed, 0=leftmost)
 * in the computational basis, leaving other n_qubits-1 qubits in renormalized
 * superposition survives outcome
 *
 * P(outcome=0) = sum over all basis states i with bit `target` = 0 of
 * |\psi[i]|^2. Samples outcome in {0,1} using caller-supplied uniform random u
 * in [0,1). Collapses \psi in place: amplitudes inconsistent with sampled
 * outcome are zeroed, survivors are rescaled by 1/\sqrt(P(outcome)) so state
 * remains normalized.
 *
 * psi->n must equal 2^n_qubits (caller-supplied n_qubits, since state vector
 * itself carries no qubit-count metadata).
 *
 * Returns sampled outcome (0 or 1), or -1 on invalid input (NULL psi, target
 * out of [0, n_qubits), or psi->n != 2^n_qubits).
 */
int qstate_measure_qubit(cvector_t *psi, int n_qubits, int target, double u);

#endif
