/*
 * Five-qubit "perfect" quantum error-correcting code (Reference: Laflamme,
 * Miquel, Paz, Zurek 1996; Bennett, DiVincenzo, Smolin, Wootters 1996)
 *
 * Encodes one logical qubit into five physical qubits and corrects an
 * arbitrary single-qubit error X, Y, or Z on any one of them. It is called
 * "perfect" because it saturates quantum Hamming bound: with n = 5 physical
 * qubits, 4 stabilizer generators and 2^n = 32 basis states, 16 syndrome
 * classes (1 no-error + 15 single-qubit Paulis) exactly fill space, leaving
 * no room for a second error and none wasted
 *
 * Stabilizers are four cyclic shifts of XZZXI over data qubits:
 *
 *     S1 = X Z Z X I
 *     S2 = I X Z Z X
 *     S3 = X I X Z Z
 *     S4 = Z X I X Z
 *
 * Code space is simultaneous +1 eigenspace of these four commuting operators
 * Measuring all four gives a 4-bit syndrome that uniquely identifies which of
 * 15 single-qubit errors (if any) occurred; syndrome-to-correction map is
 * qec5_syndrome_table in this file
 *
 * Rather than constructing code space via stabilizer formalism, two codewords
 * |0_L>, |1_L> are hard-coded as amplitude patterns over 32 computational basis
 * states. Each is a uniform superposition of 16 states with amplitudes of
 * magnitude 1/4 and signs listed in qec5_codeword[]. These signs are correct
 * Z-eigenvalue pattern for chosen stabilizers, and test suite verifies that
 * state obtained by running explicit ancilla circuit on |0_L> really is in +1
 * eigenspace of every S_k. This keeps encoder free of stabilizer machinery
 * while remaining exact
 *
 * Syndrome extraction uses standard phase-kickback circuit: prepare each
 * ancilla in |+>, apply controlled-Pauli for each non-identity letter of
 * stabilizer word, apply H, measure in Z. Ancillas q5..q8 are real qubits on
 * same state vector, so entire routine is just a fixed sequence of single- and
 * two-qubit gates plus four projective measurements
 *
 * Readout: after correction, ancillas hold measured syndrome bits in their low
 * 4 bits, so logical amplitudes are read at index (codeword << 4) | anc_val.
 * There is no reset-to-|0> step because correction is applied directly to data
 * qubits without touching ancillas and decode sums over data-qubit subspace at
 * fixed ancilla pattern
 */

#include "qec5.h"
#include "../core/complex.h"
#include "../core/vector.h"
#include "angular.h"
#include "physics/qec.h"
#include "qubits.h"
#include <math.h>
#include <stddef.h>

// Stabilizer generators S1..S4: cyclic shifts of X Z Z X I over data qubits
// 0..4
static const char *const qec5_stabilizers[4] = {"XZZXI", "IXZZX", "XIXZZ",
                                                "ZXIXZ"};

// |0_L> nonzero amplitudes: 16 of 32 five-qubit basis states, each with
// magnitude 1/4 and sign given here (index = 5-bit pattern q0q1q2q3q4,
// MSB-first)
// |1_L> = X_L|0_L> has identical signs at bitwise 5-bit complement (31-index)
// of each entry
static const struct {
  int index;
  int sign;
} qec5_codeword[16] = {
    {0, +1},  {3, -1},  {5, +1},  {6, -1},  {9, +1},  {10, +1},
    {12, -1}, {15, -1}, {17, -1}, {18, +1}, {20, +1}, {23, -1},
    {24, -1}, {27, -1}, {29, -1}, {30, -1},
};

// Syndrome -> (errored qubit, Pauli type) lookup, indexed by
// s = 8*s1 + 4*s2 + 2*s3 + s4 (s_k = measured outcome of ancilla for S_k)
// Derived from stabilizers' commutation relations with each of 15 possible
// single-qubit errors, and independently cross-checked by simulating exact
// ancilla H/controlled-Pauli/H/measure circuit below: all 16 syndromes (1 "no
// error" + 15 single-qubit errors) come out distinct, confirming code
// saturates quantum Hamming bound
static const struct {
  int qubit;
  qec_error_type_t type;
} qec5_syndrome_table[16] = {
    /* s= 0, 0000 */ {-1, 0},
    /* s= 1, 0001 */ {0, QEC_ERROR_X},
    /* s= 2, 0010 */ {2, QEC_ERROR_Z},
    /* s= 3, 0011 */ {4, QEC_ERROR_X},
    /* s= 4, 0100 */ {4, QEC_ERROR_Z},
    /* s= 5, 0101 */ {1, QEC_ERROR_Z},
    /* s= 6, 0110 */ {3, QEC_ERROR_X},
    /* s= 7, 0111 */ {4, QEC_ERROR_Y},
    /* s= 8, 1000 */ {1, QEC_ERROR_X},
    /* s= 9, 1001 */ {3, QEC_ERROR_Z},
    /* s=10, 1010 */ {0, QEC_ERROR_Z},
    /* s=11, 1011 */ {0, QEC_ERROR_Y},
    /* s=12, 1100 */ {2, QEC_ERROR_X},
    /* s=13, 1101 */ {1, QEC_ERROR_Y},
    /* s=14, 1110 */ {2, QEC_ERROR_Y},
    /* s=15, 1111 */ {3, QEC_ERROR_Y},
};

static const complex_t *pauli_gate(qec_error_type_t type) {
  switch (type) {
  case QEC_ERROR_X:
    return sigma_x;
  case QEC_ERROR_Y:
    return sigma_y;
  case QEC_ERROR_Z:
  default:
    return sigma_z;
  }
}

qec5_result_t qec5_run(complex_t alpha, complex_t beta, int error_qubit,
                       qec_error_type_t error_type, const double u[4]) {
  qec5_result_t result = {0};
  result.corrected_qubit = -1;

  double norm = sqrt(c_abs2(alpha) + c_abs2(beta));
  if (norm < 1e-300 || !u) {
    return result;
  }

  alpha = c_scale(alpha, 1.0 / norm);
  beta = c_scale(beta, 1.0 / norm);

  const int n = 9; // q0..q4 = data, q5..q8 = ancillas (one per stabilizer)

  cvector_t *psi = qstate_alloc(n);
  for (long long i = 0; i < psi->n; i++) {
    psi->data[i] = c_zero();
  }

  // Encode: direct amplitude assignment
  // Ancillas start at |0000>, so every populated basis state has its low 4 bits
  // (ancilla qubits 5..8) equal to 0, i.e. full index = data_index << 4
  for (int k = 0; k < 16; k++) {
    int idx0 = qec5_codeword[k].index;
    double sign = (double)qec5_codeword[k].sign;
    int idx1 = 31 - idx0; // |1_L>'s amplitude sits at bit-complement

    psi->data[idx0 << 4] =
        c_add(psi->data[idx0 << 4], c_scale(alpha, sign * 0.25));
    psi->data[idx1 << 4] =
        c_add(psi->data[idx1 << 4], c_scale(beta, sign * 0.25));
  }

  // Inject error
  if (error_qubit >= 0 && error_qubit < 5) {
    qstate_apply_gate1(psi, n, error_qubit, pauli_gate(error_type));
  }

  // Syndrome extraction: one ancilla-based stabilizer measurement per S_k,
  // via standard phase-kickback circuit (H, controlled-Pauli's, H, measure)
  for (int k = 0; k < 4; k++) {
    int anc = 5 + k;

    qstate_apply_gate1(psi, n, anc, hadamard_gate);
    for (int q = 0; q < 5; q++) {
      char p = qec5_stabilizers[k][q];
      if (p == 'I') {
        continue;
      }

      const complex_t *gate = (p == 'X') ? sigma_x : sigma_z;
      qstate_apply_controlled_u(psi, n, anc, q, gate);
    }

    qstate_apply_gate1(psi, n, anc, hadamard_gate);

    result.syndrome[k] = qstate_measure_qubit(psi, n, anc, u[k]);
  }

  int s = (result.syndrome[0] << 3) | (result.syndrome[1] << 2) |
          (result.syndrome[2] << 1) | result.syndrome[3];

  result.corrected_qubit = qec5_syndrome_table[s].qubit;
  result.corrected_type = qec5_syndrome_table[s].type;

  if (result.corrected_qubit >= 0) {
    qstate_apply_gate1(psi, n, result.corrected_qubit,
                       pauli_gate(result.corrected_type));
  }

  // Decode: ancillas now hold fixed measured syndrome bits (low 4 bits of every
  // surviving basis-state index); project data-qubit subspace at that fixed
  // ancilla pattern onto |0_L>/|1_L> codewords directly
  int anc_val = s;
  complex_t rec_alpha = c_zero();
  complex_t rec_beta = c_zero();

  for (int k = 0; k < 16; k++) {
    int idx0 = qec5_codeword[k].index;
    double sign = (double)qec5_codeword[k].sign;
    int idx1 = 31 - idx0;

    rec_alpha = c_add(rec_alpha,
                      c_scale(psi->data[(idx0 << 4) | anc_val], sign * 0.25));
    rec_beta =
        c_add(rec_beta, c_scale(psi->data[(idx1 << 4) | anc_val], sign * 0.25));
  }

  result.recovered_alpha = rec_alpha;
  result.recovered_beta = rec_beta;

  cvector_free(psi);

  return result;
}
