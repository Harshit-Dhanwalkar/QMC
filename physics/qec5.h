#ifndef QMC_QEC5_H
#define QMC_QEC5_H

#include "../core/complex.h"
#include "qec.h" // reuses qec_error_type_t (QEC_ERROR_X/Y/Z)

/*
 * 5-qubit "perfect" code ([[5,1,3]], References: Laflamme-Miquel-Paz-Zurek 1996
 * / Bennett-DiVincenzo-Smolin-Wootters 1996): the smallest stabilizer code that
 * corrects an arbitrary single-qubit error (X, Y, or Z, on any one of the 5
 * physical qubits).
 * "Perfect" refers to the code saturating the quantum Hamming bound: 2^4 = 16
 * distinct syndromes exactly cover the 1 "no error" outcome + 5 qubits * 3
 * Pauli types = 15 single-qubit error outcomes, with none left over and none
 * colliding.
 *
 * Stabilizer generators (cyclic shifts of X Z Z X I over the 5 data qubits):
 *   S1 = X0 Z1 Z2 X3 I4
 *   S2 = I0 X1 Z2 Z3 X4
 *   S3 = X0 I1 X2 Z3 Z4
 *   S4 = Z0 X1 I2 X3 Z4
 * Logical operators: X_L = X0 X1 X2 X3 X4, Z_L = Z0 Z1 Z2 Z3 Z4.
 *
 * NOTE: 9-qubit system: qubits 0-4 are the 5 data (physical) qubits; qubits 5-8
 * are ancillas, one dedicated per stabilizer check (S1..S4 respectively).
 * - Encoding: unlike qec_run/qec_shor_run, the logical codewords |0_L>, |1_L>
 * are prepared by direct amplitude assignment from the stabilizer projector P =
 * (1/16) * prod_k (I + S_k) applied to |00000>, rather than via an explicit
 * CNOT/H encoding circuit: each codeword is an equal-magnitude superposition of
 * 16 of the 32 five-qubit basis states (amplitude +-1/4, signs from
 * qec5_codeword_table below), and |1_L> = X_L|0_L> which flips every data-qubit
 * bit, so its nonzero support is exactly the bitwise 5-bit-complement of
 * |0_L>'s support with identical signs. The physically parts of the protocol
 * error injection, syndrome extraction (ancilla H / controlled-Pauli / H /
 * measure phase-kickback circuit, one round per stabilizer), and correction -
 * are gate circuits
 * - Syndrome extraction (per stabilizer S_k, ancilla a_k initialized to |0>):
 *    H(a_k); for each data qubit i with S_k having a non-identity Pauli P at
 *    i, apply controlled-P(a_k -> i); H(a_k); measure a_k in Z basis.
 *    Outcome 0 => +1 eigenvalue (stabilizer satisfied on that check), 1 => -1.
 * The resulting 4-bit syndrome (s1,s2,s3,s4) is looked up in
 * qec5_syndrome_table to identify errored qubit and Pauli type, or (-1, -1) for
 * no error.
 * - Correction: unlike 9-qubit Shor code (which infers X and Z corrections
 * independently, so a physical Y error gets corrected as Z*X, recovering the
 * state up to an unavoidable global phase of i), the 5-qubit code's syndrome
 * identifies the exact error type (X, Y, or Z) directly, so recovery applies
 * that exact Pauli and is exact with no phase ambiguity for any of the 15
 * possible single-qubit errors.
 * - Decoding: since encoding bypassed an explicit circuit, decoding likewise
 * reads out the logical amplitudes by direct inner product against the
 * |0_L>/|1_L> codewords (<0_L|psi>, <1_L|psi>) restricted to the data-qubit
 * subspace at whatever fixed ancilla bit pattern the post-measurement state now
 * occupies, rather than by an inverse circuit.
 */

typedef struct {
  int syndrome[4]; // measured ancilla bits s1..s4 (one per stabilizer S1..S4)
  int corrected_qubit; // data qubit (0..4) identified as errored, or
                       // -1 if syndrome indicated no error
  qec_error_type_t corrected_type; // Pauli type applied as correction
  complex_t recovered_alpha;       // decoded logical qubit after correction
  complex_t recovered_beta; // recovered_alpha*|0> + recovered_beta*|1>, should
                            // exactly equal original (\alpha, \beta) regardless
                            // of which single-qubit error (if any) was injected
} qec5_result_t;

/*
 * Encode \alpha*|0> + \beta*|1> (normalized internally) into the 5-qubit code,
 * on `error_qubit`, run 4-round ancilla syndrome-extraction circuit, apply the
 * indicated correction, and decode.
 *
 * NOTE: u[4] supplies caller randomness for the 4 ancilla measurements (one per
 * stabilizer), same caller-supplied-randomness convention as
 * qstate_measure_qubit; as in qec_run, every measurement is deterministic
 * (probability exactly 0 or 1) given the injected error, so any u[k] in [0,1)
 * yields the same correct result.
 */
qec5_result_t qec5_run(complex_t alpha, complex_t beta, int error_qubit,
                       qec_error_type_t error_type, const double u[4]);

#endif
