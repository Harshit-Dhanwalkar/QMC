/*
Quantum teleportation, superdense coding, CHSH/Bell inequality.
*/

#include "quantum_info.h"
#include "../core/complex.h"
#include "../core/vector.h"
#include "angular.h"
#include "qubits.h"
#include <math.h>
#include <stdlib.h>

qi_teleport_result_t qi_teleport(complex_t alpha, complex_t beta, double u1,
                                 double u2) {
  qi_teleport_result_t result = {0};

  /* 3 qubits: q0 = \alpha|0> + \beta|1> (to be teleported),
   *  q1 = Alice's half of a fresh Bell pair,
   *  q2 = Bob's half.
   *  qstate_alloc starts at |000>;
   *  set q0's amplitudes directly (normalize defensively).
   */
  double norm = sqrt(c_abs2(alpha) + c_abs2(beta));
  if (norm < 1e-300) {
    return result;
  }

  alpha = c_scale(alpha, 1.0 / norm);
  beta = c_scale(beta, 1.0 / norm);

  cvector_t *psi = qstate_alloc(3);
  psi->data[0] = alpha; // |000> : q0=0
  psi->data[4] = beta;  // |100> : q0=1 (qubit 0 = leftmost/MSB)

  // Prepare Bell pair on q1,q2
  qstate_apply_gate1(psi, 3, 1, hadamard_gate);
  qstate_apply_cnot(psi, 3, 1, 2);

  // Alice's Bell-basis measurement circuit on q0,q1
  qstate_apply_cnot(psi, 3, 0, 1);
  qstate_apply_gate1(psi, 3, 0, hadamard_gate);

  int m1 = qstate_measure_qubit(psi, 3, 0, u1);
  int m2 = qstate_measure_qubit(psi, 3, 1, u2);

  // After both measurements, psi is a product state |m1>|m2>|bob>; read Bob's
  // (q2) amplitudes directly
  int base = (m1 << 2) | (m2 << 1);
  complex_t bob_alpha = psi->data[base + 0];
  complex_t bob_beta = psi->data[base + 1];

  // Correction: X if m2==1, then Z if m1==1
  if (m2 == 1) {
    complex_t tmp = bob_alpha;
    bob_alpha = bob_beta;
    bob_beta = tmp;
  }
  if (m1 == 1) {
    bob_beta = c_scale(bob_beta, -1.0);
  }

  cvector_free(psi);

  result.m1 = m1;
  result.m2 = m2;
  result.bob_alpha = bob_alpha;
  result.bob_beta = bob_beta;

  return result;
}

qi_superdense_result_t qi_superdense(int bit1, int bit2) {
  qi_superdense_result_t result = {0};

  cvector_t *psi = qstate_alloc(2);

  // Shared Bell pair |Phi+>
  qstate_apply_gate1(psi, 2, 0, hadamard_gate);
  qstate_apply_cnot(psi, 2, 0, 1);

  // Alice encodes on her qubit (q0): 00->I, 01->X, 10->Z, 11->X then Z
  if (bit1 == 0 && bit2 == 1) {
    qstate_apply_gate1(psi, 2, 0, sigma_x);
  } else if (bit1 == 1 && bit2 == 0) {
    qstate_apply_gate1(psi, 2, 0, sigma_z);
  } else if (bit1 == 1 && bit2 == 1) {
    qstate_apply_gate1(psi, 2, 0, sigma_x);
    qstate_apply_gate1(psi, 2, 0, sigma_z);
  }

  // Bob decodes
  qstate_apply_cnot(psi, 2, 0, 1);
  qstate_apply_gate1(psi, 2, 0, hadamard_gate);

  // Deterministic by construction: exactly one basis state has amplitude ~1,
  // rest ~0
  int best = 0;
  double best_p = 0.0;
  for (int i = 0; i < 4; i++) {
    double p = c_abs2(psi->data[i]);
    if (p > best_p) {
      best_p = p;
      best = i;
    }
  }

  cvector_free(psi);

  result.decoded_bit1 = (best >> 1) & 1;
  result.decoded_bit2 = best & 1;

  return result;
}

// A(\theta) = \cos(\theta) * sigma_z + \sin(\theta) * \sigma_x, as a 2x2 gate
// array
static void spin_measurement_operator(double theta, complex_t out[4]) {
  double c = cos(theta), s = sin(theta);
  for (int k = 0; k < 4; k++) {
    out[k] = c_add(c_scale(sigma_z[k], c), c_scale(sigma_x[k], s));
  }
}

double qi_chsh_correlator(double theta, double phi) {
  complex_t A[4], B[4];
  spin_measurement_operator(theta, A);
  spin_measurement_operator(phi, B);

  // Bell state |\Phi+> = (|00> + |11>) / \sqrt2
  cvector_t *psi = cvector_alloc(4);
  psi->data[0] = c_real(M_SQRT1_2);
  psi->data[1] = c_zero();
  psi->data[2] = c_zero();
  psi->data[3] = c_real(M_SQRT1_2);

  // (A tensor B) * \psi, via explicit Kronecker-product structure:
  // (A tensor // B)[2i + k][2j + l] = A[i][j] * B[k][l]
  cvector_t *op_psi = cvector_alloc(4);
  for (int i = 0; i < 4; i++) {
    op_psi->data[i] = c_zero();
  }

  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 2; k++) {
      complex_t sum = c_zero();
      for (int j = 0; j < 2; j++) {
        for (int l = 0; l < 2; l++) {
          complex_t coeff = c_mul(A[i * 2 + j], B[k * 2 + l]);

          sum = c_add(sum, c_mul(coeff, psi->data[j * 2 + l]));
        }
      }

      op_psi->data[i * 2 + k] = sum;
    }
  }

  double e = cvector_expect(psi, op_psi);

  cvector_free(psi);
  cvector_free(op_psi);

  return e;
}

double qi_chsh_S(double a, double a_prime, double b, double b_prime) {
  return qi_chsh_correlator(a, b) - qi_chsh_correlator(a, b_prime) +
         qi_chsh_correlator(a_prime, b) + qi_chsh_correlator(a_prime, b_prime);
}
