/*
3-qubit bit-flip and phase-flip QEC codes, ancilla-based syndrome extraction.
*/

#include "qec.h"
#include "../core/complex.h"
#include "../core/vector.h"
#include "angular.h"
#include "qubits.h"
#include <math.h>

qec_result_t qec_run(qec_code_t code, complex_t alpha, complex_t beta,
                     int error_qubit, double u3, double u4) {
  qec_result_t result = {0};

  result.corrected_qubit = -1;

  double norm = sqrt(c_abs2(alpha) + c_abs2(beta));
  if (norm < 1e-300) {
    return result;
  }

  alpha = c_scale(alpha, 1.0 / norm);
  beta = c_scale(beta, 1.0 / norm);

  const int n = 5; // q0,q1,q2 = data, q3,q4 = ancilla

  cvector_t *psi = qstate_alloc(n);
  // |00000>: index 0. |11100> (q0=q1=q2=1, q3=q4=0): bit weights 16+8+4=28
  // (qubit k's bit is at position n-1-)
  int idx000 = 0;
  int idx111 = 28;
  psi->data[idx000] = alpha;
  psi->data[idx111] = beta;

  const complex_t *error_gate = (code == QEC_BITFLIP) ? sigma_x : sigma_z;

  if (code == QEC_PHASEFLIP) {
    qstate_apply_gate1(psi, n, 0, hadamard_gate);
    qstate_apply_gate1(psi, n, 1, hadamard_gate);
    qstate_apply_gate1(psi, n, 2, hadamard_gate);
  }

  if (error_qubit >= 0 && error_qubit < 3) {
    qstate_apply_gate1(psi, n, error_qubit, error_gate);
  }

  if (code == QEC_PHASEFLIP) {
    qstate_apply_gate1(psi, n, 0, hadamard_gate);
    qstate_apply_gate1(psi, n, 1, hadamard_gate);
    qstate_apply_gate1(psi, n, 2, hadamard_gate);
  }

  // Syndrome extraction: Z0Z1 and Z1Z2 parity onto two ancillas
  qstate_apply_cnot(psi, n, 0, 3);
  qstate_apply_cnot(psi, n, 1, 3);
  qstate_apply_cnot(psi, n, 1, 4);
  qstate_apply_cnot(psi, n, 2, 4);

  if (code == QEC_PHASEFLIP) {
    qstate_apply_gate1(psi, n, 0, hadamard_gate);
    qstate_apply_gate1(psi, n, 1, hadamard_gate);
    qstate_apply_gate1(psi, n, 2, hadamard_gate);
  }

  int s3 = qstate_measure_qubit(psi, n, 3, u3);
  int s4 = qstate_measure_qubit(psi, n, 4, u4);

  int corrected_qubit;
  if (s3 == 0 && s4 == 0) {
    corrected_qubit = -1;
  } else if (s3 == 1 && s4 == 0) {
    corrected_qubit = 0;
  } else if (s3 == 1 && s4 == 1) {
    corrected_qubit = 1;
  } else { // s3==0 && s4==1
    corrected_qubit = 2;
  }

  if (corrected_qubit >= 0) {
    qstate_apply_gate1(psi, n, corrected_qubit, error_gate);
  }

  if (code == QEC_PHASEFLIP) {
    qstate_apply_gate1(psi, n, 0, hadamard_gate);
    qstate_apply_gate1(psi, n, 1, hadamard_gate);
    qstate_apply_gate1(psi, n, 2, hadamard_gate);
  }

  // Decoded logical amplitudes: ancillas now hold measured syndrome (s3,s4),
  // not |00>, so readout index must include their actual post-measurement
  // values
  int idx000_actual = (s3 << 1) | s4;
  int idx111_actual = 28 | (s3 << 1) | s4;

  result.syndrome_s3 = s3;
  result.syndrome_s4 = s4;
  result.corrected_qubit = corrected_qubit;
  result.recovered_alpha = psi->data[idx000_actual];
  result.recovered_beta = psi->data[idx111_actual];

  cvector_free(psi);

  return result;
}

// Bit-flip syndrome extraction + correction for one block {a,b,c}, reusing
// ancillas anc1/anc2 (which must be |0> on entry and are reset back to |0> on
// exit).
// Returns corrected qubit index, or -1 if syndrome indicated no error. u[0],
// u[1] supply two ancilla measurements.
static int shor_bitflip_correct_block(cvector_t *psi, int n, int a, int b,
                                      int c, int anc1, int anc2,
                                      const double u[2]) {
  qstate_apply_cnot(psi, n, a, anc1);
  qstate_apply_cnot(psi, n, b, anc1);
  qstate_apply_cnot(psi, n, b, anc2);
  qstate_apply_cnot(psi, n, c, anc2);

  int s1 = qstate_measure_qubit(psi, n, anc1, u[0]);
  int s2 = qstate_measure_qubit(psi, n, anc2, u[1]);

  int corrected;
  if (s1 == 0 && s2 == 0) {
    corrected = -1;
  } else if (s1 == 1 && s2 == 0) {
    corrected = a;
  } else if (s1 == 1 && s2 == 1) {
    corrected = b;
  } else {
    corrected = c;
  }

  if (corrected >= 0) {
    qstate_apply_gate1(psi, n, corrected, sigma_x);
  }

  // Reset ancillas back to |0> for reuse by the next block/check.
  if (s1 == 1) {
    qstate_apply_gate1(psi, n, anc1, sigma_x);
  }
  if (s2 == 1) {
    qstate_apply_gate1(psi, n, anc2, sigma_x);
  }

  return corrected;
}

qec_shor_result_t qec_shor_run(complex_t alpha, complex_t beta, int error_qubit,
                               qec_error_type_t error_type, const double u[8]) {
  qec_shor_result_t result = {0};
  result.block_corrected[0] = result.block_corrected[1] =
      result.block_corrected[2] = -1;
  result.phase_block_corrected = -1;

  double norm = sqrt(c_abs2(alpha) + c_abs2(beta));
  if (norm < 1e-300 || !u) {
    return result;
  }

  alpha = c_scale(alpha, 1.0 / norm);
  beta = c_scale(beta, 1.0 / norm);

  const int n = 11; // q0..q8 = data (3 blocks of 3), q9,q10 = ancillas
  const int blocks[3][3] = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};

  cvector_t *psi = qstate_alloc(n);
  // |0>_L -> all-zero index 0; |1>_L -> qubit 0 set, everything else 0:
  // bit weight 1 << (n-1-0) = 1 << 10 = 1024.
  psi->data[0] = alpha;
  psi->data[1 << (n - 1)] = beta;

  // Encode
  qstate_apply_cnot(psi, n, 0, 3);
  qstate_apply_cnot(psi, n, 0, 6);
  qstate_apply_gate1(psi, n, 0, hadamard_gate);
  qstate_apply_gate1(psi, n, 3, hadamard_gate);
  qstate_apply_gate1(psi, n, 6, hadamard_gate);
  qstate_apply_cnot(psi, n, 0, 1);
  qstate_apply_cnot(psi, n, 0, 2);
  qstate_apply_cnot(psi, n, 3, 4);
  qstate_apply_cnot(psi, n, 3, 5);
  qstate_apply_cnot(psi, n, 6, 7);
  qstate_apply_cnot(psi, n, 6, 8);

  // Inject error
  if (error_qubit >= 0 && error_qubit < 9) {
    const complex_t *gate = (error_type == QEC_ERROR_X)   ? sigma_x
                            : (error_type == QEC_ERROR_Y) ? sigma_y
                                                          : sigma_z;
    qstate_apply_gate1(psi, n, error_qubit, gate);
  }

  // Bit-flip correction, one block at a time, reusing ancillas 9,10
  for (int b = 0; b < 3; b++) {
    result.block_corrected[b] = shor_bitflip_correct_block(
        psi, n, blocks[b][0], blocks[b][1], blocks[b][2], 9, 10, &u[2 * b]);
  }

  // Phase-flip correction across blocks: stabilizers X0..X5, X3..X8 (all 6
  // qubits of each pair of blocks being compared, not just their leaders)
  for (int q = 0; q < 9; q++) {
    qstate_apply_gate1(psi, n, q, hadamard_gate);
  }
  for (int q = 0; q < 6; q++) {
    qstate_apply_cnot(psi, n, q, 9);
  }
  for (int q = 3; q < 9; q++) {
    qstate_apply_cnot(psi, n, q, 10);
  }
  for (int q = 0; q < 9; q++) {
    qstate_apply_gate1(psi, n, q, hadamard_gate);
  }

  int ps1 = qstate_measure_qubit(psi, n, 9, u[6]);
  int ps2 = qstate_measure_qubit(psi, n, 10, u[7]);

  int phase_block;
  if (ps1 == 0 && ps2 == 0) {
    phase_block = -1;
  } else if (ps1 == 1 && ps2 == 0) {
    phase_block = 0;
  } else if (ps1 == 1 && ps2 == 1) {
    phase_block = 1;
  } else {
    phase_block = 2;
  }

  if (phase_block >= 0) {
    qstate_apply_gate1(psi, n, blocks[phase_block][0], sigma_z);
  }
  if (ps1 == 1) {
    qstate_apply_gate1(psi, n, 9, sigma_x);
  }
  if (ps2 == 1) {
    qstate_apply_gate1(psi, n, 10, sigma_x);
  }

  result.phase_block_corrected = phase_block;

  // Decode: exact inverse of encoding circuit
  qstate_apply_cnot(psi, n, 6, 8);
  qstate_apply_cnot(psi, n, 6, 7);
  qstate_apply_cnot(psi, n, 3, 5);
  qstate_apply_cnot(psi, n, 3, 4);
  qstate_apply_cnot(psi, n, 0, 2);
  qstate_apply_cnot(psi, n, 0, 1);
  qstate_apply_gate1(psi, n, 0, hadamard_gate);
  qstate_apply_gate1(psi, n, 3, hadamard_gate);
  qstate_apply_gate1(psi, n, 6, hadamard_gate);
  qstate_apply_cnot(psi, n, 0, 6);
  qstate_apply_cnot(psi, n, 0, 3);

  // Qubits 1..8 are back to |0>; qubits 3,6 are back to |0> too. Ancillas 9,10
  // hold their last measured values (already reset to |0> above), so logical
  // amplitudes sit at index 0 (|0>_L) and 1<<(n-1) (|1>_L), exactly as at
  // encoding time
  result.recovered_alpha = psi->data[0];
  result.recovered_beta = psi->data[1 << (n - 1)];

  cvector_free(psi);

  return result;
}
