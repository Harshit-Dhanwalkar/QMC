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

#endif
