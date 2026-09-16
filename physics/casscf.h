#ifndef QMC_CASSCF_H
#define QMC_CASSCF_H

#include "../core/matrix.h"
#include "molecular_hf.h"
#include "molecular_integrals.h"

/*
 * NOTE: CASSCF (Complete Active Space Self-Consistent Field): simultaneously
 * optimizes the CI coefficients and the molecular orbitals themselves, unlike
 * CISD/FCI (physics/cisd.h, physics/fci.h) which diagonalize a fixed-orbital
 * Hamiltonian. Partitions the n_basis spatial orbitals into n_frozen (always
 * doubly occupied, excluded from the CI expansion), n_active (a full CI within
 * this active space), and the remainder (n_basis - n_frozen - n_active, always
 * empty "virtual" orbitals) - exactly same partition used by CASSCF in every
 * quantum chemistry package.
 *
 * Algorithm:
 *   1. Active-space CI: build the frozen-core-effective Hamiltonian over just
 *      the n_frozen+n_active orbitals, diagonalize in the n_electrons_active
 *      sector, and construct the 1- and 2-particle reduced density matrices
 *      (RDMs) from the resulting CI vector via direct bit-manipulation
 *      fermionic operators
 *   2. Extend the active-space RDMs to the FULL orbital space (core
 *      orbitals contribute a fixed closed-shell block; active contributes
 *      CI density; virtual is zero), using MCSCF density partition formulas
 *   3. Orbital gradient: rather than an analytically-derived generalized Fock
 *      matrix, dE/dkappa_pq for every independent orbital-rotation generator
 *      kappa_pq is computed by Central finite difference, holding the RDMs
 *      from step 1 fixed and only recomputing the (no re-diagonalization is
 *      needed) energy expression: E = E_nuc + \sum h_pq D_pq + 1/2 \sum (pq|rs)
 *      d_pqrs under a small orbital rotation. The diagonal of the Hessian is
 *      estimated the  same way (3-point finite difference), giving an
 *      approximate Newton step size per generator; the denominator is floored
 *      in absolute value (not sign-preserved) to guarantee a valid descent
 *      direction even where the true curvature is tiny or slightly negative
 *   4. Backtracking line search on the actual (re-diagononalized) active-space
 *      CI energy, then rotate C <- C @ exp(kappa) using a  real matrix
 *      exponential (scaling-and-squaring + Taylor series).
 *   5. Repeat from step 1 with the rotated orbitals until the gradient norm
 *      falls below conv_tol.
 */

typedef struct {
  double total_energy; /* electronic + nuclear repulsion */
  int n_spatial;
  int n_frozen;
  int n_active;
  int n_electrons_active;
  cmatrix_t *C;     /* optimized n_spatial x n_spatial MO coefficients */
  int iterations;   /* number of orbital-rotation steps taken */
  int converged;    /* 1 if final gradient norm < conv_tol, else 0 */
  double grad_norm; /* final orbital-gradient norm */
} casscf_result_t;

/*
 * basis/n_basis/mol: as in molecular_hf.h/molecular_integrals.h.
 * C_initial: starting MO coefficients (n_basis x n_basis), e.g. from
 *   molecular_rhf's result->C -- CASSCF is a local optimization and, like
 *   any MCSCF method, needs a reasonable starting guess.
 * n_frozen, n_active, n_electrons_active: active-space partition (see
 *   above). Requires 0 <= n_frozen, n_active >= 1,
 *   n_frozen + n_active <= n_basis, and
 *   0 <= n_electrons_active <= 2*n_active.
 * conv_tol: convergence threshold on the orbital-gradient norm (a
 *   reasonable default is 1e-6).
 * max_iter: maximum orbital-rotation iterations.
 *
 * Returns NULL on invalid input or allocation failure. Check
 * result->converged; a non-converged result still holds the best C/energy
 * found within max_iter iterations.
 */
casscf_result_t *casscf_run(basis_function_t **basis, int n_basis,
                            const molecule_t *mol, const cmatrix_t *C_initial,
                            int n_frozen, int n_active, int n_electrons_active,
                            double conv_tol, int max_iter);

void casscf_result_free(casscf_result_t *res);

#endif
