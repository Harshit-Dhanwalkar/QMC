/*
Distance-3 rotated surface code ([[9,1,3]]): smallest instance topological
surface code family
*/

#include "qec_surface17.h"
#include "../core/complex.h"
#include "../core/vector.h"
#include "angular.h"
#include "physics/qec.h"
#include "qubits.h"
#include <math.h>
#include <stddef.h>

typedef struct {
  char type; // 'X' or 'Z'
  int n_support;
  int support[4]; // data qubits 0..8 this stabilizer acts on
} qec_surface17_stabilizer_t;

// 8 stabilizer generators S0..S7
static const qec_surface17_stabilizer_t qec_surface17_stabilizers[8] = {
    {'X', 2, {3, 6}},       {'Z', 2, {0, 1}},       {'X', 4, {0, 1, 3, 4}},
    {'Z', 4, {3, 4, 6, 7}}, {'Z', 4, {1, 2, 4, 5}}, {'X', 4, {4, 5, 7, 8}},
    {'Z', 2, {7, 8}},       {'X', 2, {2, 5}},
};

// |0_L> nonzero amplitudes: 16 of 512 nine-qubit basis states, each with
// amplitude exactly +1/4
// |1_L> = X_L|0_L> (X_L = X0 X1 X2) has identical amplitudes at each index with
// data qubits 0,1,2 toggled (MSB-first bit mask 0x1C0 for a 9-qubit index)
static const int qec_surface17_codeword[16] = {
    0, 27, 36, 63, 72, 83, 108, 119, 399, 404, 427, 432, 455, 476, 483, 504,
};

#define QEC_SURFACE17_XL_MASK 0x1C0 // bits for data qubits 0,1,2 (MSB-first)

// Syndrome -> canonical (qubit, Pauli type) decode table, indexed by s =
// \sum_k(syndrome[k] << (7-k)). Derived by simulating exact ancilla
// syndrome-extraction circuit for every one of 28 possible "no error or
// single-qubit error" cases; some pairs of distinct single-qubit errors alias
// to same syndrome (expected for a distance-3, non-"perfect" code - differ by a
// weight-2 stabilizer element)
static const struct {
  int syndrome;
  int qubit;
  qec_error_type_t type;
} qec_surface17_decode_table[24] = {
    {0, -1, 0},
    {1, 2, QEC_ERROR_Z},
    {2, 8, QEC_ERROR_X},
    {4, 7, QEC_ERROR_Z},
    {5, 5, QEC_ERROR_Z},
    {6, 8, QEC_ERROR_Y},
    {8, 2, QEC_ERROR_X},
    {9, 2, QEC_ERROR_Y},
    {13, 5, QEC_ERROR_Y},
    {16, 3, QEC_ERROR_X},
    {18, 7, QEC_ERROR_X},
    {22, 7, QEC_ERROR_Y},
    {24, 4, QEC_ERROR_X},
    {32, 0, QEC_ERROR_Z},
    {36, 4, QEC_ERROR_Z},
    {60, 4, QEC_ERROR_Y},
    {64, 0, QEC_ERROR_X},
    {72, 1, QEC_ERROR_X},
    {96, 0, QEC_ERROR_Y},
    {104, 1, QEC_ERROR_Y},
    {128, 6, QEC_ERROR_Z},
    {144, 6, QEC_ERROR_Y},
    {160, 3, QEC_ERROR_Z},
    {176, 3, QEC_ERROR_Y},
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

// Linear search over 24-entry decode table; a syndrome outside this set cannot
// arise from a weight<=1 error (only kind this function ever injects), so this
// always finds a match in normal use. Falls back to "no correction" defensively
// if it somehow doesn't
static void decode_syndrome(int s, int *qubit, qec_error_type_t *type) {
  for (int i = 0; i < 24; i++) {
    if (qec_surface17_decode_table[i].syndrome == s) {
      *qubit = qec_surface17_decode_table[i].qubit;
      *type = qec_surface17_decode_table[i].type;

      return;
    }
  }

  *qubit = -1;
  *type = QEC_ERROR_X;
}

qec_surface17_result_t qec_surface17_run(complex_t alpha, complex_t beta,
                                         int error_qubit,
                                         qec_error_type_t error_type,
                                         const double u[8]) {
  qec_surface17_result_t result = {0};
  result.corrected_qubit = -1;

  double norm = sqrt(c_abs2(alpha) + c_abs2(beta));
  if (norm < 1e-300 || !u) {
    return result;
  }

  alpha = c_scale(alpha, 1.0 / norm);
  beta = c_scale(beta, 1.0 / norm);

  const int n = 17; // q0..q8 = data, q9..q16 = ancillas (one per stabilizer)

  cvector_t *psi = qstate_alloc(n);
  for (long long i = 0; i < psi->n; i++) {
    psi->data[i] = c_zero();
  }

  // Encode: direct amplitude assignment. Ancillas start at |00000000>, so every
  // populated basis state has its low 8 bits (ancilla qubits 9..16) equal to 0,
  // i.e. full index = data_index << 8
  for (int k = 0; k < 16; k++) {
    int idx0 = qec_surface17_codeword[k];
    int idx1 = idx0 ^ QEC_SURFACE17_XL_MASK;

    psi->data[idx0 << 8] = c_add(psi->data[idx0 << 8], c_scale(alpha, 0.25));
    psi->data[idx1 << 8] = c_add(psi->data[idx1 << 8], c_scale(beta, 0.25));
  }

  // Inject error
  if (error_qubit >= 0 && error_qubit < 9) {
    qstate_apply_gate1(psi, n, error_qubit, pauli_gate(error_type));
  }

  // Syndrome extraction: one ancilla-based stabilizer measurement per S_k.
  for (int k = 0; k < 8; k++) {
    int anc = 9 + k;
    const qec_surface17_stabilizer_t *stab = &qec_surface17_stabilizers[k];
    const complex_t *gate = (stab->type == 'X') ? sigma_x : sigma_z;

    qstate_apply_gate1(psi, n, anc, hadamard_gate);
    for (int i = 0; i < stab->n_support; i++) {
      qstate_apply_controlled_u(psi, n, anc, stab->support[i], gate);
    }

    qstate_apply_gate1(psi, n, anc, hadamard_gate);

    result.syndrome[k] = qstate_measure_qubit(psi, n, anc, u[k]);
  }

  int s = 0;
  for (int k = 0; k < 8; k++) {
    s |= result.syndrome[k] << (7 - k);
  }

  decode_syndrome(s, &result.corrected_qubit, &result.corrected_type);

  if (result.corrected_qubit >= 0) {
    qstate_apply_gate1(psi, n, result.corrected_qubit,
                       pauli_gate(result.corrected_type));
  }

  // Decode: ancillas now hold fixed measured syndrome bits (low 8 bits of every
  // surviving basis-state index); data-qubit subspace at that fixed ancilla
  // pattern onto |0_L>/|1_L> codewords directly
  int anc_val = s;
  complex_t rec_alpha = c_zero();
  complex_t rec_beta = c_zero();

  for (int k = 0; k < 16; k++) {
    int idx0 = qec_surface17_codeword[k];
    int idx1 = idx0 ^ QEC_SURFACE17_XL_MASK;

    rec_alpha =
        c_add(rec_alpha, c_scale(psi->data[(idx0 << 8) | anc_val], 0.25));
    rec_beta = c_add(rec_beta, c_scale(psi->data[(idx1 << 8) | anc_val], 0.25));
  }

  result.recovered_alpha = rec_alpha;
  result.recovered_beta = rec_beta;

  cvector_free(psi);

  return result;
}
