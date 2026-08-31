#ifndef QMC_FCI_H
#define QMC_FCI_H

#include "../core/matrix.h"

/*
 * Full Configuration Interaction (FCI): exact diagonalization of molecular
 * electronic Hamiltonian within a fixed-particle-number sector of Fock space,
 * given MO-basis one- and two-electron integrals.
 *
 * NOTE: This module instead explicitly restricts to requested electron-number
 * sector before diagonalizing:
 *   - Correctness independent of "does the right answer happen to also be
 *     the global minimum": works even for open-shell, high-symmetry, or
 *     dissociating cases where a different particle-number sector could in
 *     principle have a lower eigenvalue than the requested one.
 *   - Efficiency: diagonalizing a C(2*n_spatial, n_electrons)-dimensional
 *     sector matrix (e.g. 70x70 for H4/STO-3G's 4-electron-in-8-spin-
 *     orbital sector) instead of the full 2^(2*n_spatial)-dimensional
 *     space (256x256 for the same system) is both faster and uses less
 *     memory for the actual diagonalization step.
 *
 * HACK: Basis-state convention (shared with second_quant.h): Fock-space index i
 * (0 <= i < 2^(2 * n_spatial)) is a bitstring occupation number over
 * 2*n_spatial spin orbitals; electron count = popcount(i). This module's sector
 * restriction is exactly "keep only i with popcount(i) == n_electrons".
 */

typedef struct {
  double ground_energy;    /* lowest eigenvalue in the requested sector,
                            * electronic + nuclear repulsion (matches
                            * second_quant_build_molecular_hamiltonian's
                            * convention: total molecular energy, directly
                            * comparable to RHF/CCSD/etc. total energies) */
  double *eigenvalues;     /* full sector spectrum, ascending, length dim */
  cmatrix_t *eigenvectors; /* dim x dim, column k is the eigenvector for
                            * eigenvalues[k], expressed in sector basis (i.e.
                            * eigenvectors->data row index r corresponds to
                            * Fock-space state basis_states[r], not r itself) */
  int *basis_states;       /* length dim: the Fock-space index (bitstring) each
                            * sector-basis row/column corresponds to, so callers can
                            * map back to full Fock-space occupation patterns if
                            * needed (e.g. to compute properties, or identify dominant
                            * configuration) */
  int dim;                 /* C(2*n_spatial, n_electrons) */
  int n_spatial;
  int n_electrons;
} fci_result_t;

/*
 * Full CI within the n_electrons-electron sector, using all n_spatial
 * spatial (2*n_spatial spin) orbitals as active space.
 *
 * h_mo, eri_mo, nuclear_repulsion: identical convention to
 * second_quant_build_molecular_hamiltonian (h_mo is n_spatial^2 one-electron
 * integrals, eri_mo is n_spatial^4 two-electron integrals in chemist
 * notation, nuclear_repulsion is added as a constant energy shift).
 *
 * Returns NULL if n_spatial<1, h_mo/eri_mo is NULL, or n_electrons is not
 * in [0, 2*n_spatial].
 */
fci_result_t *fci_solve(int n_spatial, const double *h_mo, const double *eri_mo,
                        double nuclear_repulsion, int n_electrons);

/*
 * Frozen-core FCI: freezes first n_frozen spatial orbitals as doubly occupied
 * (excluded from active space and correlation treatment entirely), then runs
 * FCI with n_electrons_active electrons in remaining (n_spatial-n_frozen)
 * active spatial orbitals. Uses
 * second_quant_build_molecular_hamiltonian_frozen_core internally; In
 * second_quant.h for exact frozen-core energy-shift and effective-Hamiltonian
 * formulas.
 *
 * Returns NULL under the same conditions as
 * second_quant_build_molecular_hamiltonian_frozen_core, or if
 * n_electrons_active is not in [0, 2*(n_spatial-n_frozen)].
 */
fci_result_t *fci_solve_frozen_core(int n_spatial, int n_frozen,
                                    const double *h_mo, const double *eri_mo,
                                    double nuclear_repulsion,
                                    int n_electrons_active);

void fci_result_free(fci_result_t *res);

#endif
