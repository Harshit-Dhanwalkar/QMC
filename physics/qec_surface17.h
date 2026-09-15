#ifndef QMC_QEC_SURFACE17_H
#define QMC_QEC_SURFACE17_H

#include "../core/complex.h"
#include "qec.h" // reuses qec_error_type_t (QEC_ERROR_X/Y/Z)

/*
 * Distance-3 rotated surface code: a [[9,1,3]] CSS stabilizer code (9 data
 * qubits, 1 logical qubit, distance 3), theoretical smallest instance of
 * topological surface code family (Reference: Fowler et al., "Surface codes:
 * Towards practical large-scale quantum computation", 2012; same
 * 9-data/8-ancilla layout underlies Google's "Surface-17" experiments)
 * Unlike qec5 (which corrects every single-qubit error exactly, saturating
 * quantum Hamming bound), a distance-3 surface code only guarantees correcting
 * errors of weight <= 1; some pairs of distinct single-qubit errors alias to
 * same syndrome (they differ by a weight-2 stabilizer element), which is
 * expected and harmless; correcting with either labeled qubit's Pauli exactly
 * restores logical state.
 *
 * - Lattice (data qubits 0-8, row-major on a 3x3 grid):
 *   0  1  2
 *   3  4  5
 *   6  7  8
 * - Stabilizers (8 generators: 4 X-type "plaquette", 4 Z-type "vertex"),
 *   derived from rotated-surface-code construction; plaquette centers on dual
 *   lattice, checkerboard-colored, weight-4 in bulk, weight-2 on boundary with
 *   only "matching color" boundary plaquettes kept (X-type weight-2 on
 *   left/right edges, Z-type weight-2 on top/bottom edges), verified to commute
 *   pairwise, be linearly independent (GF(2) rank 8), and yield a 2-dimensional
 *   (1-logical-qubit) code space:
 *   S0 = X3 X6            (weight 2, left edge)
 *   S1 = Z0 Z1            (weight 2, top edge)
 *   S2 = X0 X1 X3 X4      (weight 4, bulk)
 *   S3 = Z3 Z4 Z6 Z7      (weight 4, bulk)
 *   S4 = Z1 Z2 Z4 Z5      (weight 4, bulk)
 *   S5 = X4 X5 X7 X8      (weight 4, bulk)
 *   S6 = Z7 Z8            (weight 2, bottom edge)
 *   S7 = X2 X5            (weight 2, right edge)
 *   Logical operators (brute-force confirmed minimum weight = 3 = code
 *     distance): X_L = X0 X1 X2 (top row), Z_L = Z0 Z3 Z6 (left column)
 *     17-qubit system: qubits 0-8 are 9 data qubits; qubits 9-16 are ancillas,
 *     one dedicated per stabilizer check (S0..S7 respectively)
 * - Encoding: as with qec5, |0_L>/|1_L> are prepared by Direct amplitude *
 *   assignment from stabilizer projector P = (1/256) * prod_k (I + S_k) applied
 *   to |000000000>, rather than an explicit encoding circuit. Because only 4
 *   X-type stabilizers move |000000000> (Z-type ones fix it exactly), codeword
 *   collapses to just 16 (not 256) nonzero basis-state amplitudes, each +1/4.
 *   |1_L> = X_L|0_L> flips data qubits {0,1,2}, so its support is exactly
 *   qec_surface17_codeword's indices with bits 0,1,2 (MSB-first) toggled, at
 *   same sign (+1)
 * - Syndrome extraction: identical ancilla-based protocol to qec5 (H,
 *   controlled-Pauli's matching each stabilizer's type/support, H, measure),
 *   one round per stabilizer, using real controlled-X/controlled-Z gates:
 *   explicit gate circuit
 * - Decoding: syndrome is looked up in qec_surface17_decode_table to get
 *   canonical (qubit, Pauli type) correction (or "no error"); that Pauli is
 *   applied, and logical amplitudes are read out by direct inner product
 *   against |0_L>/|1_L> codewords
 */

typedef struct {
  int syndrome[8]; // measured ancilla bits s0..s7 (one per stabilizer S0..S7)
  int corrected_qubit; // data qubit (0..8) identified as errored, or -1
                       // if syndrome indicated no error
  qec_error_type_t corrected_type; // Pauli type applied as correction
                                   // (meaningful only if corrected_qubit>=0)
  complex_t recovered_alpha;       // decoded logical qubit after correction
  complex_t recovered_beta;        // recovered_alpha*|0> + recovered_beta*|1>
} qec_surface17_result_t;

/*
 * Encode alpha*|0> + beta*|1> (normalized internally) into distance-3 surface
 * code, optionally inject a single-qubit error of type `error_type` on
 * `error_qubit` (0..8; pass error_qubit=-1 for no error), run 8-round ancilla
 * syndrome-extraction circuit, apply indicated correction, and decode
 *
 * u[8] supplies caller randomness for 8 ancilla measurements, every measurement
 * is deterministic given a weight<=1 injected error, so any u[k] in [0,1)
 * yields same correct result
 */
qec_surface17_result_t qec_surface17_run(complex_t alpha, complex_t beta,
                                         int error_qubit,
                                         qec_error_type_t error_type,
                                         const double u[8]);

#endif
