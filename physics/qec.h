#ifndef QMC_QEC_H
#define QMC_QEC_H

#include "../core/complex.h"

/*
 * 3-qubit quantum error correction codes (bit-flip and phase-flip),
 * with syndrome extraction via two ancilla qubits.
 * NOTE: The syndrome (two ancilla measurement outcomes) reveals which data
 * qubit (if any) suffered an error without collapsing encoded logical qubit's
 * \alpha/\beta superposition, which is entire point of ancilla-based approach:
 * a direct Z-basis measurement of data qubits would destroy very superposition
 * the code is protecting.
 *
 * 5-qubit system: qubits 0,1,2 are data (logical) qubits, qubits 3,4 are
 * ancillas used only for syndrome extraction.
 *
 * BIT-FLIP CODE (protects against X errors):
 *   Encode  : \alpha|0> + \beta|1> -> \alpha|000> + \beta|111>
 *             (CNOT(0,1),CNOT(0,2))
 *   Syndrome: CNOT(0,3), CNOT(1,3), CNOT(1,4), CNOT(2,4),
 * then measure ancillas 3,4 in Z basis.
 * NOTE: This extracts Z0Z1 and Z1Z2 parity stabilizers onto ancillas: measuring
 * an ancilla doesn't collapse which computational-basis-pair superposition (000
 * vs 111, possibly with one bit flipped) data is in, only WHICH pair.
 * Syndrome table: (0,0)->no error, (1,0)->qubit 0, (1,1)->qubit 1,
 *                 (0,1)->qubit 2.
 *
 * Phase-flip code (protects against Z errors): identical circuit, conjugated by
 * Hadamard on each data qubit before encoding/after decoding and around the
 * syndrome extraction (H * X * H = Z, so bit-flip code's X-error protection
 * becomes Z-error protection in Hadamard-conjugated basis.
 */

typedef enum { QEC_BITFLIP = 0, QEC_PHASEFLIP = 1 } qec_code_t;

typedef struct {
  int syndrome_s3, syndrome_s4; // measured ancilla bits
  int corrected_qubit; // data qubit (0,1,2) identified as errored, or -1 if
                       // syndrome indicated no error
  complex_t recovered_alpha; // decoded logical qubit after correction
  complex_t recovered_beta;  // recovered_alpha*|0> + recovered_beta*|1>, should
                            // exactly equal original (\alpha, \beta) regardless
                            // of which single data qubit (if any) was hit by
                            // injected error
} qec_result_t;

/*
 * Encode \alpha * |0> + \beta * |1> (normalized internally), optionally inject
 * single-qubit error on `error_qubit` (0, 1, or 2; pass -1 for no error : an X
 * error for QEC_BITFLIP, a Z error for QEC_PHASEFLIP), run syndrome-extraction
 * circuit, apply indicated correction, and decode.
 *
 * NOTE: u3, u4 are passed to ancilla measurements using same
 * caller-supplied-randomness convention as qstate_measure_qubit; in this
 * noise-free demonstration, ancilla outcome probabilities are always exactly 0
 * or 1 given injected error (no genuine randomness is involved, unlike e.g. a
 * Bell measurement), so any valid u in [0,1) gives same, correct result. */
qec_result_t qec_run(qec_code_t code, complex_t alpha, complex_t beta,
                     int error_qubit, double u3, double u4);

/*
 * 9-qubit Shor code: corrects an ARBITRARY single-qubit error (X, Y, or Z on
 * any one of the 9 physical qubits) by concatenating the bit-flip and
 * phase-flip codes above - 3 blocks of the 3-qubit bit-flip code, with the  3
 * block "values" additionally protected by an outer phase-flip structure.
 *
 * NOTE: 11-qubit system: qubits 0-8 are the 9 physical data qubits, grouped
 * into 3 blocks {0,1,2}, {3,4,5}, {6,7,8}; qubits 9,10 are ancillas, reused
 * across all 4 syndrome extractions (3 bit-flip + 1 phase-flip), matching now a
 * real device would measure-and-reset ancillas between rounds rather than
 * dedicating fresh ones to each check.
 * - Encode: CNOT(0,3), CNOT(0,6) [copy to block leaders], H(0), H(3), H(6)
 * [rotate leaders to phase-protecting basis], then CNOT(0,1), CNOT(0,2),
 * CNOT(3,4), CNOT(3,5), CNOT(6,7), CNOT(6,8) [bit-flip-encode each block].
 * - Bit-flip correction (one call per block, using qec_run's own bit-flip
 * syndrome table): for block {a,b,c}, CNOT(a,anc1), CNOT(b,anc1), CNOT(b,anc2),
 * CNOT(c,anc2), measure both ancillas, apply X to the indicated qubit (if any),
 * then reset both ancillas back to |0> (apply X to any ancilla that measured
 * 1).
 * - Phase-flip correction (comparing block PARITIES, not just leaders - needs
 * stabilizers X0X1X2X3X4X5 and X3X4X5X6X7X8, all 6 qubits of the two blocks
 * being compared, not just their leaders, since a block's logical value is
 * carried jointly by all 3 of its qubits once bit-flip-encoded): H on all 9
 * data qubits, CNOT from block-1's 3 qubits into anc1, CNOT from block-2's 3
 * qubits into anc1 (X0..X5 parity), CNOT from block-2's 3 qubits into anc2,
 * CNOT from block-3's 3 qubits into anc2 (X3..X8 parity), H on all 9 data
 * qubits again, measure both ancillas, apply Z to a representative qubit of the
 * indicated block (if any).
 * - Decode: the exact inverse of the encoding circuit (same gates, reverse
 * order - CNOT and H are self-inverse, so this is well-defined), then read off
 * qubit 0's amplitudes.
 *
 * WARN: Y errors: correcting a Y error necessarily means the bit-flip syndrome
 * fires (treating it as an X error) AND the phase-flip syndrome fires (treating
 * it as a Z error) independently, so the applied correction is Z*X, not Y
 * itself. Since Y = iXZ, this recovers the correct logical state up to an
 * unavoidable, physically unobservable global phase of i (Z*X*Y = i*Identity) -
 * this is standard, expected behavior for independent Pauli-frame correction,
 * not an approximation and recovered_alpha/recovered_beta will reflect that
 * exact global phase for Y errors specifically (X and Z errors, and no error,
 * recover the exact original alpha/beta with no phase ambiguity at all).
 *
 * error_qubit in [0,8] selects which physical qubit gets the injected error
 * (error_type selects X, Y, or Z); pass error_qubit=-1 for no error. u[8]
 * supplies caller randomness for the 8 ancilla measurements (3 bit-flip blocks
 * x 2 ancillas each, + 1 phase-flip check x 2 ancillas) as with qec_run, every
 * measurement here is deterministic (probability exactly 0 or 1) given the
 * injected error, so any u in [0,1) per slot gives the same correct result.
 */
typedef enum {
  QEC_ERROR_X = 0,
  QEC_ERROR_Y = 1,
  QEC_ERROR_Z = 2
} qec_error_type_t;

typedef struct {
  int block_corrected[3];    // which qubit (if any, else -1) was bit-flip
                             // corrected in each of the 3 blocks
  int phase_block_corrected; // which block (0, 1, or 2), if any (-1
                             // otherwise), was phase-corrected
  complex_t recovered_alpha;
  complex_t recovered_beta;
} qec_shor_result_t;

qec_shor_result_t qec_shor_run(complex_t alpha, complex_t beta, int error_qubit,
                               qec_error_type_t error_type, const double u[8]);

#endif
