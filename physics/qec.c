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
