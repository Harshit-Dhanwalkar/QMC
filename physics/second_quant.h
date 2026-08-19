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

/*
 * Build the full second-quantized molecular electronic Hamiltonian (dense, 2^(2
 * * n_spatial) x 2^(2 * n_spatial)) via same direct bit-manipulation
 * Jordan-Wigner construction as second_quant_build_hopping_hamiltonian (not via
 * jw_creation_operator/jw_annihilation_operator tensor products), on 2 *
 * n_spatial spin-orbitals using the interleaved convention: spin-orbital index
 * 2p = spatial MO p, spin-up (\alpha); 2p+1 = spatial MO p, spin-down (\beta).
 *
 *   H = nuclear_repulsion * I
 *       + \sum_{pq} h_pq a_p^\dagger a_q
 *       + (1/2) * \sum_{pqrs} <pq|rs> a_p^\dagger a_q^\dagger a_s a_r
 *
 * h_mo             : n_spatial x n_spatial one-electron (kinetic + nuclear
 *                    attraction) integrals in the MO basis, row-major
 *                    (h_mo[p * n_spatial + q]).
 * eri_mo           : n_spatial^4 two-electron integrals in MO basis, CHEMIST
 *                    notation (pq|rs) with same flat indexing as
 *                    molecular_integrals.h's MOLINT_ERI macro (i.e. pass
 *                    exactly one get from AO->MO-transforming
 *                    molecular_eri_tensor's output). This function internally
 *                    converts to physicist <pq|rs> = chemist (pr|qs) ordering
 *                    sum above needs.
 * nuclear_repulsion: added as a constant shift to every diagonal element, so
 *                    the returned matrix's eigenvalues are directly total
 *                    molecular energies (electronic + nuclear), not just
 *                    electronic ones.
 *
 * Returns NULL if n_spatial < 1 or h_mo/eri_mo is NULL.
 */
cmatrix_t *second_quant_build_molecular_hamiltonian(int n_spatial,
                                                    const double *h_mo,
                                                    const double *eri_mo,
                                                    double nuclear_repulsion);

/*
 * Same as second_quant_build_molecular_hamiltonian, but freezes the first
 * n_frozen spatial orbitals (indices 0..n_frozen-1 in h_mo/eri_mo) as
 * doubly-occupied and excludes them from the active space entirely :
 * frozen-core approximation (e.g. freezing a tightly-bound, chemically inert
 * atomic core orbital like Li's 1s in LiH). Builds the Hamiltonian on only the
 * remaining n_spatial-n_frozen active spatial orbitals (2*(n_spatial-n_frozen)
 * spin-orbitals), with the frozen orbitals' contribution folded into (a) a
 * constant energy shift (E_core, added on top of nuclear_repulsion) and (b) an
 * effective one-electron Hamiltonian for the active orbitals that includes the
 * frozen orbitals' mean-field (Hartree + exchange) effect:
 *
 *   E_core = nuclear_repulsion + \sum_c 2 h_cc
 *            + \sum_{c,d} [2 (cc|dd) - (cd|dc)]     (c,d range over frozen)
 *   h_eff_pq = h_pq + \sum_c [2 (pq|cc) - (pc|cq)]  (c ranges over frozen, p,q
 *                                                    over active)
 *
 * with the active-space two-electron integrals otherwise unchanged (just
 * restricted to active-active-active-active indices).
 *
 * Returns NULL if n_spatial < 1, n_frozen < 0, n_frozen >= n_spatial, or
 * h_mo/eri_mo is NULL.
 */
cmatrix_t *second_quant_build_molecular_hamiltonian_frozen_core(
    int n_spatial, int n_frozen, const double *h_mo, const double *eri_mo,
    double nuclear_repulsion);

#endif
