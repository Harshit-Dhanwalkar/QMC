/*
 * Deutsch-Jozsa and Grover's search algorithms
 *
 * 2 theoritical query-complexity algorithms implemented on top of dense
 * state-vector simulator in qubits.c. Both assume an N-qubit register stored as
 * a length-2^N complex vector psi->data, with qubit index 0 as most-significant
 * bit of state index
 *
 *  1. deutsch_jozsa()
 *     One-shot determination of whether a boolean oracle f: {0,1}^n -> {0,1} is
 *     constant or balanced, using standard phase-kickback circuit: ancilla in
 *     |1> gets Hadamarded into |->, inputs are put into a uniform
 *     superposition, oracle is applied once, then inputs are Hadamarded back
 *     Result p_all_zero is probability of measuring all input qubits as 0; for
 *     a constant oracle this is 1, for a balanced oracle it is 0 (in ideal
 *     limit). is_constant is derived from p_all_zero > 0.5
 *     Oracle is specified via dj_oracle_type_t:
 *       DJ_CONSTANT_0  identity (f = 0)
 *       DJ_CONSTANT_1  ancilla flip (f = 1)
 *       DJ_BALANCED    CNOTs from a caller-chosen subset of input qubits
 *                      onto ancilla (f = XOR of those bits)
 *  2. grover_search()
 *     Amplitude-amplification search over 2^n_qubits basis states for a
 *     caller-specified target index. Each iteration applies oracle (phase flip
 *     on target amplitude) followed by diffusion operator 2|s><s| - I expressed
 *     as H^n (2|0><0| - I) H^n. Default iteration count is round((\pi/4)
 *     \sqrt(N)) via lround. Returns full probability distribution over N basis
 *     states along with number of iterations actually performed
 */

#include "quantum_algorithms.h"
#include "../core/complex.h"
#include "../core/vector.h"
#include "angular.h"
#include "qubits.h"
#include <math.h>
#include <stdlib.h>

dj_result_t deutsch_jozsa(int n_input, dj_oracle_type_t oracle_type,
                          const int *parity_qubits, int n_parity_qubits) {
  dj_result_t result = {-1.0, 0};

  if (n_input < 1) {
    return result;
  }
  if (oracle_type == DJ_BALANCED && (!parity_qubits || n_parity_qubits < 1)) {
    return result;
  }

  int n = n_input + 1; // + ancilla, last qubit (index n_input)
  int ancilla = n_input;

  cvector_t *psi = qstate_alloc(n); // (or n_qubits)
  if (!psi) {
    return result; // (deutsch_jozsa: result already zeroed)
  }

  // Prepare ancilla in |1>: qstate_alloc starts at |0...0>, so flip ancilla,
  // then Hadamard every qubit (inputs -> uniform superposition, ancilla -> |->)
  qstate_apply_gate1(psi, n, ancilla, sigma_x);
  for (int q = 0; q < n; q++) {
    qstate_apply_gate1(psi, n, q, hadamard_gate);
  }

  // Oracle (phase kickback via ancilla in |->)
  if (oracle_type == DJ_CONSTANT_1) {
    qstate_apply_gate1(psi, n, ancilla, sigma_x);
  } else if (oracle_type == DJ_BALANCED) {
    for (int k = 0; k < n_parity_qubits; k++) {
      qstate_apply_cnot(psi, n, parity_qubits[k], ancilla);
    }
  }
  // DJ_CONSTANT_0: oracle is identity, nothing to apply

  for (int q = 0; q < n_input; q++) {
    qstate_apply_gate1(psi, n, q, hadamard_gate);
  }

  // P(all input qubits = 0): with ancilla as last (lowest-order) qubit, "input
  // qubits all zero" means state-vector index's high bits (everything except
  // ancilla bit) are zero, i.e. index < 2
  double p0 = c_abs2(psi->data[0]) + c_abs2(psi->data[1]);

  cvector_free(psi);

  result.p_all_zero = p0;
  result.is_constant = (p0 > 0.5) ? 1 : 0;

  return result;
}

grover_result_t grover_search(int n_qubits, int target, int n_iterations) {
  grover_result_t result = {0};

  if (n_qubits < 1) {
    return result;
  }

  int N = 1 << n_qubits;
  if (target < 0 || target >= N) {
    return result;
  }

  if (n_iterations < 0) {
    n_iterations = (int)lround((M_PI / 4.0) * sqrt((double)N));
  }

  cvector_t *psi = qstate_alloc(n_qubits);
  for (int q = 0; q < n_qubits; q++) {
    qstate_apply_gate1(psi, n_qubits, q, hadamard_gate);
  }

  for (int it = 0; it < n_iterations; it++) {
    // Oracle: flip phase of target amplitude
    psi->data[target] = c_scale(psi->data[target], -1.0);

    // NOTE: Diffusion operator 2|s><s| - I = H^n (2|0><0|-I) H^n: Hadamard
    // every qubit, negate every amplitude except index 0 (this is exactly
    // 2|0><0|-I: diagonal +1 at index 0, -1 elsewhere), Hadamard every qubit
    // again
    for (int q = 0; q < n_qubits; q++) {
      qstate_apply_gate1(psi, n_qubits, q, hadamard_gate);
    }
    for (int i = 1; i < N; i++) {
      psi->data[i] = c_scale(psi->data[i], -1.0);
    }
    for (int q = 0; q < n_qubits; q++) {
      qstate_apply_gate1(psi, n_qubits, q, hadamard_gate);
    }
  }

  double *probabilities = malloc((size_t)N * sizeof *probabilities);
  if (!probabilities) {
    cvector_free(psi);

    return result;
  }

  for (int i = 0; i < N; i++) {
    probabilities[i] = c_abs2(psi->data[i]);
  }

  result.n_qubits = n_qubits;
  result.target = target;
  result.n_iterations = n_iterations;
  result.p_target = probabilities[target];
  result.probabilities = probabilities;

  cvector_free(psi);

  return result;
}
