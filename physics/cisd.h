#ifndef QMC_CISD_H
#define QMC_CISD_H

#include "../core/matrix.h"

/*
 * Configuration Interaction Singles and Doubles (CISD): variational
 * diagonalization of the molecular electronic Hamiltonian restricted to
 * reference (canonical HF) determinant plus all single and double excitations
 * out of it, rather than the full N-electron sector (FCI). Natural extension of
 * fci.h/second_quant.h's infrastructure: reuses the same full Fock-space
 * Hamiltonian builder and just restricts to a smaller basis-state subset before
 * diagonalizing (Reference: Szabo & Ostlund Ch. 4.3).
 *
 * Reference determinant: spin-orbitals 0..n_electrons-1 occupied (lowest
 * interleaved indices, see second_quant.h's 2p=alpha_p/2p+1=beta_p convention),
 * i.e. the canonical HF ground-state occupation IF h_mo/eri_mo are passed in
 * increasing-orbital-energy (canonical HF) order
 * NOTE: second_quant.h's Fock-space index has spin-orbital 0 as the MSB
 * (leftmost) bit (qstate_*'s convention), so this reference is  bit pattern
 * with the TOP n_electrons bits of the n_modes-bit index set.
 *
 * Basis-state selection (bit-manipulation criterion, shared with fci.h's
 * Fock-space index convention: bitstring i, bit j set means spin-orbital j
 * occupied): a determinant D (popcount(D) == n_electrons) is included iff
 * popcount(D XOR reference) is 0 (the reference itself), 2 (a single
 * excitation: exactly one occupied->virtual swap), or 4 (a double excitation:
 * exactly two occupied->virtual swaps). This dimension formula:
 *  1 + n_occ * n_virt + C(n_occ,2) * C(n_virt,2)
 * with n_occ=n_electrons, n_virt=2 * n_spatial-n_electrons.
 *
 * CISD is variational (an eigenvalue of a Hermitian matrix restricted to a
 * subspace containing reference), so ground_energy always satisfies
 *  E_reference >= E_CISD >= E_FCI
 * for same system/active space - a sanity check against fci.h/hartree_fock.h
 * results.
 */

typedef struct {
  double ground_energy;      /* lowest eigenvalue in the CISD subspace,
                              * electronic + nuclear repulsion (directly
                              * comparable to RHF/FCI/CCSD total energies) */
  double reference_energy;   /* <D0|H|D0>, the reference determinant's own
                              * diagonal Hamiltonian matrix element (should
                              * match the input RHF total energy if the
                              * integrals came from a converged RHF/UHF
                              * calculation and the reference is its ground
                              * state) */
  double correlation_energy; /* ground_energy - reference_energy (<= 0) */
  double *eigenvalues;       /* full CISD-subspace spectrum, ascending,
                              * length dim */
  cmatrix_t *eigenvectors;   /* dim x dim, column k is the eigenvector for
                              * eigenvalues[k], expressed in CISD-subspace
                              * basis (row index r <-> basis_states[r]) */
  int *basis_states;         /* length dim: Fock-space index (bitstring) each
                              * CISD-subspace row/column corresponds to */
  int dim;                   /* 1 + n_occ * n_virt + C(n_occ,2) * C(n_virt,2) */
  int n_spatial;
  int n_electrons;
} cisd_result_t;

/*
 * CISD within the n_electrons-electron sector, using all n_spatial spatial
 * (2 * n_spatial spin) orbitals as active space.
 *
 * h_mo, eri_mo, nuclear_repulsion: identical convention to fci_solve /
 * second_quant_build_molecular_hamiltonian.
 *
 * Returns NULL if n_spatial<1, h_mo/eri_mo is NULL, or n_electrons is not
 * in [0, 2*n_spatial].
 */
cisd_result_t *cisd_solve(int n_spatial, const double *h_mo,
                          const double *eri_mo, double nuclear_repulsion,
                          int n_electrons);

/*
 * Frozen-core CISD: identical frozen-core treatment to fci_solve_frozen_core
 * (freezes the first n_frozen spatial orbitals as doubly occupied, excluded
 * entirely from the active space and singles/doubles excitation manifold), then
 * runs CISD with n_electrons_active electrons in remaining (n_spatial-n_frozen)
 * active spatial orbitals.
 *
 * Returns NULL under the same conditions as fci_solve_frozen_core.
 */
cisd_result_t *cisd_solve_frozen_core(int n_spatial, int n_frozen,
                                      const double *h_mo, const double *eri_mo,
                                      double nuclear_repulsion,
                                      int n_electrons_active);

void cisd_result_free(cisd_result_t *res);

#endif
