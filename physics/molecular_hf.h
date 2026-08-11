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
