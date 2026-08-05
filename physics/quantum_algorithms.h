#ifndef QMC_QUANTUM_ALGORITHMS_H
#define QMC_QUANTUM_ALGORITHMS_H

#include "../core/vector.h"

/*
 * 2 canonical quantum algorithms :
 * - Deutsch-Jozsa's constant/balanced cases give exactly P=1/P=0,
 * - Grover's search matched closed-form success-probability formula :
 *      \sin^2((2k + 1) * \theta), \sin(\theta)=1 / \sqrt(N),
 * to ~1e-6 across several (n_qubits, target) combinations.
 *
 * NOTE:
 * - Oracles are implemented as direct phase operations on state vector (e.g.
 * negating a single target amplitude for Grover's oracle) rather than
 * decomposed into elementary multi-controlled gates
 * - TFIM Hamiltonian builder uses constructing intended operator directly
 * rather than via a full gate decomposition. A real quantum computer would need
 * decomposed circuit; a state-vector simulator does not.
 */

typedef enum {
  DJ_CONSTANT_0 = 0, // f(x) = 0 for all x
  DJ_CONSTANT_1 = 1, // f(x) = 1 for all x
  DJ_BALANCED = 2    // f(x) = XOR of input bits listed in parity_qubits (1/2
                     // 2^n_input inputs give 0, 1/2 give 1)
} dj_oracle_type_t;

typedef struct {
  double p_all_zero; // probability of measuring all n_input input qubits as 0
  int is_constant; // algorithm's conclusion: 1 if p_all_zero > 1/2 (constant),
                   // 0 otherwise (balanced) with a noise-free oracle always
                   // exactly 1.0 or 0.0
} dj_result_t;

/*
 * Deutsch-Jozsa: determines whether an oracle f: {0,1}^n_input -> {0,1} is
 * constant or balanced with a single oracle query (vs. up to 2^(n_input - 1) +
 * 1 classical queries in worst case).
 *
 * For DJ_BALANCED, parity_qubits (length n_parity_qubits, each entry a valid
 * input-qubit index in [0, n_input)) selects which input bits oracle XORs
 * together; f(x) = XOR of those bits is balanced for any non-empty subset.
 *
 * Returns a zeroed dj_result_t (p_all_zero = -1.0) on invalid input (n_input <
 * 1, or a DJ_BALANCED oracle with no parity qubits given).
 */
dj_result_t deutsch_jozsa(int n_input, dj_oracle_type_t oracle_type,
                          const int *parity_qubits, int n_parity_qubits);

typedef struct {
  int n_qubits;
  int target;
  int n_iterations;
  double p_target; // probability of measuring target state after n_iterations
                   // Grover iterations
  double *probabilities; // full probability distribution, length 2^n_qubits
} grover_result_t;

/*
 * Grover's search: amplifies amplitude of a single marked `target` state (0 <=
 * target < 2^n_qubits) among an initially uniform superposition.
 *
 * Returns a zeroed grover_result_t (probabilities == NULL) on invalid input
 * (n_qubits < 1, target out of range).
 */
grover_result_t grover_search(int n_qubits, int target, int n_iterations);

#endif
