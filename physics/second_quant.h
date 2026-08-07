#ifndef QMC_SECOND_QUANT_H
#define QMC_SECOND_QUANT_H

#include "../core/matrix.h"

/*
 * Second quantization and the Jordan-Wigner (JW) transformation: mapping
 * fermionic creation/annihilation operators to qubit operators, so fermionic
 * many-body Hamiltonians can be built and diagonalized instead of needing a
 * dedicated fermionic-Fock-space solver.
 *
 * NOTE: Convention: for n_modes fermionic modes, 2^n_modes-dimensional Fock
 * space is identified with n_modes-qubit computational basis (mode j occupied
 * <-> qubit j is |1>, matching qubit-0-is-leftmost/MSB convention).
 * The Jordan-Wigner mapping is :
 *   a_j^\dagger = Z_0 Z_1 ... Z_{j-1} (x) \sigma^+_j (x) I_{j+1} ... I_{n-1}
 *   a_j         = Z_0 Z_1 ... Z_{j-1} (x) \sigma^-_j (x) I_{j+1} ... I_{n-1}
 *
 * Where
 *  \sigma^+ = |1><0| (creates a particle: |0>->|1>),
 *  \sigma^- = |0><1| (annihilates: |1>->|0>)
 * and leading Z-string on lower-indexed modes supplies fermionic
 * anticommutation sign (trick that makes qubit tensor-product operators satisfy
 * fermionic algebra despite qubits themselves being bosonic degrees of
 * freedom).
 */

/*
 * Jordan-Wigner creation/annihilation operator for fermionic mode `mode`
 * (0-indexed) out of n_modes total modes, as a dense 2^n_modes x 2^n_modes
 * cmatrix_t. Built via explicit Z-string x \sigma^+/- tensor-product
 * construction.
 *
 * Returns NULL if mode is out of [0, n_modes) or n_modes < 1.
 */
cmatrix_t *jw_creation_operator(int mode, int n_modes);
cmatrix_t *jw_annihilation_operator(int mode, int n_modes);

/*
 * Builds dense (2^n_modes x 2^n_modes) Hamiltonian of a spinless fermionic
 * tight-binding chain with nearest-neighbor interaction (a toy "spinless
 * Hubbard" model), directly via bit manipulation and explicit fermionic sign
 * tracking (Not via jw_creation_operator/jw_annihilation_operator or any tensor
 * products. Hamiltonian buildersg:
 *   H = \sum_i \epsilon_i n_i -
 *       t * \sum_i (a_i^\dagger a_{i+1} + a_{i+1}^\dagger a_i) +
 *       U * \sum_i n_i n_{i+1}
 *
 * Where
 *  \epsilon: length n_modes, on-site energies
 *  t       : hopping amplitude.
 *  U       : nearest-neighbor density-density interaction strength (U=0
 *            recovers a purely non-interacting tight-binding chain).
 *
 * Returns NULL if n_modes < 1 or epsilon is NULL.
 */
cmatrix_t *second_quant_build_hopping_hamiltonian(int n_modes,
                                                  const double *epsilon,
                                                  double t, double U);

#endif
