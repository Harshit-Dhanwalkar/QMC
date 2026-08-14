#ifndef QMC_MOLECULAR_HF_H
#define QMC_MOLECULAR_HF_H

#include "../core/matrix.h"
#include "molecular_integrals.h"

/*
 * General N-basis-function closed-shell restricted Hartree-Fock (RHF) and
 * AO->MO integral transform. (built on top of general GTO engine).
 */

typedef struct {
  int n_basis;
  int n_electrons;
  double total_energy; /* electronic + nuclear repulsion */
  double electronic_energy;
  cmatrix_t *C;             /* n_basis x n_basis MO coefficient matrix
                             * (columns are MOs, ascending orbital energy) */
  double *orbital_energies; /* length n_basis */
  int iterations;
  int converged;
} molecular_hf_result_t;

/*
 * Run closed-shell RHF for a general basis/molecule. n_electrons must be even
 * basis and mol as built by molecular_integrals.c (molint_basis_sto3g_h,
 * basis_function_alloc, molecule_alloc, ...). Uses symmetric (Lowdin)
 * orthogonalization via overlap matrix's own eigendecomposition, then
 * Roothaan-equation SCF iteration with simple density-based convergence.
 *
 * Returns NULL on invalid input (n_basis<=0, odd/nonpositive n_electrons,
 * n_electrons > 2*n_basis) or allocation failure.
 */
molecular_hf_result_t *molecular_rhf(basis_function_t **basis, int n_basis,
                                     const molecule_t *mol, int n_electrons,
                                     double tol, int max_iter);

void molecular_hf_result_free(molecular_hf_result_t *res);

/*
 * General N-basis-function open-shell Unrestricted Hartree-Fock (UHF):
 * separate (\alpha / \beta) density matrices and MO coefficients, additive
 * extension of molecular_rhf. Unlocks radicals, open-shell cations/anions, and
 * triplet (or higher-multiplicity) states that RHF's paired-electron assumption
 * cannot represent.
 *
 * Fock matrices (Refrence: UHF equations, Szabo & Ostlund Ch. 3.8):
 *   F^a_ij = Hcore_ij + J[D^a+D^b]_ij - K[D^a]_ij
 *   F^b_ij = Hcore_ij + J[D^a+D^b]_ij - K[D^b]_ij
 * Where
 *  D^s_ij = sum over occupied spin-s MOs of C^s_ik C^s_jk (no factor of 2 -
 * single-particle density per spin channel, unlike RHF's paired D). Both spin
 * channels are diagonalized against the same Lowdin-orthogonalized basis
 * (shared S^{-1/2}), each self-consistently, sharing only Coulomb (J) term
 * built from the total density.
 */
typedef struct {
  int n_basis;
  int n_alpha, n_beta; /* n_alpha >= n_beta by convention, but not enforced */
  double total_energy; /* electronic + nuclear repulsion */
  double electronic_energy;
  cmatrix_t *C_alpha, *C_beta; /* n_basis x n_basis MO coefficient matrices */
  double *orbital_energies_alpha; /* length n_basis */
  double *orbital_energies_beta;  /* length n_basis */
  double spin_squared;            /* <S^2> diagnostic: exactly Sz(Sz+1) for an
                                   * uncontaminated pure spin state (Sz=(n_a-n_b)/2);
                                   * values above that indicate spin contamination from
                                   * higher-multiplicity determinants mixing in */
  int iterations;
  int converged;
} molecular_uhf_result_t;

/*
 * Run UHF for a general basis/molecule. n_alpha and n_beta may differ
 * (open-shell); n_alpha,n_beta >= 0, n_alpha+n_beta >= 1,
 * max(n_alpha,n_beta) <= n_basis. tol/max_iter as in molecular_rhf.
 *
 * Returns NULL on invalid input or allocation failure.
 */
molecular_uhf_result_t *molecular_uhf(basis_function_t **basis, int n_basis,
                                      const molecule_t *mol, int n_alpha,
                                      int n_beta, double tol, int max_iter);

void molecular_uhf_result_free(molecular_uhf_result_t *res);

/*
 * NOTE: AO -> MO transform for both one-electron core Hamiltonian and
 * two-electron integral tensor, via RHF (or any other) MO coefficient matrix C.
 * h_mo must be preallocated n_basis*n_basis.
 * eri_mo must be preallocated n_basis^4 (both row-major / MOLINT_ERI
 * flat-indexed).
 *
 * WARN: since it is only ever used here on same small active spaces
 * second_quant_build_molecular_hamiltonian is already restricted to. Nnot
 * intended to scale to large basis sets without optimising this).
 */
void molecular_ao_to_mo(const cmatrix_t *h_ao, const double *eri_ao,
                        const cmatrix_t *C, int n_basis, double *h_mo,
                        double *eri_mo);

#endif
