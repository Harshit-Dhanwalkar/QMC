#ifndef QMC_TDDFT_H
#define QMC_TDDFT_H

#include "molecular_dft.h"
#include "molecular_integrals.h"

/*
 * Linear-response TDDFT within the Tamm-Dancoff Approximation (TDA), for a
 * converged closed-shell Kohn-Sham LDA ground state from molecular_dft.h.
 * Natural extension of that module's SCF infrastructure: reuses its converged
 * MO coefficients/orbital energies, AO integrals, and Becke grid
 *
 * NOTE: THEORY (References: Casida 1995; see also Dreuw & Head-Gordon, Chem.
 * Rev. 105, 4009 (2005) Sec. 2 for the TDA specialization): for a closed-shell
 * reference, singlet vertical excitation energies are eigenvalues of the (real,
 * symmetric) coupling matrix, indexed by occupied-virtual pairs (i,a):
 *
 *   A_{ia,jb} = \delta_ij \delta_ab (eps_a - eps_i) + 2 * (ia|jb)
 *               + 2 * \int \phi_i(r) \phi_a(r) f_xc(r) \phi_j(r) \phi_b(r) dr
 *
 * Where:
 *  (ia|jb) is the chemist-notation two-electron mo integral, eps are the
 * converged ks orbital energies, and f_xc(r) = d(v_xc)/dn |_{n(r)} is the
 * (local, since lda) exchange-correlation kernel, evaluated at the ground-state
 * density n(r). this formula (including both factor-of-2's, appropriate for a
 * spin-restricted singlet response with a pure, non-hybrid functional)
 *
 * f_xc(r) is computed by central finite difference of lda_xc_potential(n)
 * (slater exchange + pz81 correlation, exactly what molecular_ks_lda's SCF
 * already uses)
 *
 * WARN: Scope: tda only (not the full non-hermitian casida equations), LDA only
 * (no PBE/GGA kernel; a GGA kernel needs density-gradient response terms not
 * implemented), singlet excitations only, and excitation energies only (no
 * oscillator strengths / transition dipoles)
 */

typedef struct {
  double *excitation_energies; /* ascending, length dim, Hartree */
  int dim;                     /* n_occ * n_virt */
  int n_occ;
  int n_virt;
  int n_basis;
} tddft_result_t;

/*
 * ks must be a converged (ks->converged != 0) molecular_ks_lda result for
 * the same basis/mol; grid must be the same grid (or an equivalent one)
 * used to converge it.
 *
 * Returns NULL if ks is not converged, n_basis<1, n_electrons is odd or
 * exceeds 2*n_basis, the resulting occupied-virtual space is empty
 * (n_occ==0 or n_virt==0), or any allocation fails.
 */
tddft_result_t *tddft_tda_lda(basis_function_t **basis, int n_basis,
                              const molecule_t *mol,
                              const molecular_dft_result_t *ks,
                              const molecular_grid_t *grid);

void tddft_result_free(tddft_result_t *res);

#endif
