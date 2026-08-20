#ifndef QMC_UCC_H
#define QMC_UCC_H

#include "../core/matrix.h"
#include "../core/vector.h"
#include "ccsd.h"

/*
 * Unitary Coupled Cluster (UCC) ansatz, converged T1/T2 amplitudes as its
 * parameter source.
 *
 * |\psi(\theta)> = \exp(\kappa(\theta)) |reference>,
 * \kappa anti-Hermitian:
 *   \kappa = \sum_k \theta_s[k] * (a_a^+ a_i - a_i^+ a_a)
 *          + \sum_k \theta_d[k] * (a_a^+ a_b^+ a_j a_i - a_i^+ a_j^+ a_b a_a)
 *   i.e. \kappa = [singles] + [doubles]
 *
 * Unlike CCSD's non-unitary \exp(T) (T alone, not T-T^\dagger), this
 * exponential is manifestly unitary by construction (\kappa anti-Hermitian =>
 * \exp(\kappa) unitary), making it a legitimate variational ansatz for a
 * quantum computer.
 *
 * NOTE: dense 2^n_modes x 2^n_modes generator matrix and its diagonalization -
 * fine for the small molecules (H2, H4, LiH), not intended to scale past that
 * without a Pauli-string/Trotterized implementation.
 */

typedef struct {
  int i, a; /* occupied -> virtual spin-orbital indices */
} ucc_single_t;

typedef struct {
  int i, j, a, b; /* occupied pair -> virtual pair spin-orbital indices,
                     excitation (i,j) -> (a,b) */
} ucc_double_t;

/*
 * Build dense anti-Hermitian generator \kappa(\theta) (2^n_modes x 2^n_modes)
 * from explicit singles/doubles excitation lists and their \theta parameters.
 * Either list may be NULL if its count is 0.
 *
 * Returns NULL on invalid input (n_modes<1, or any index in singles/doubles out
 * of [0,n_modes)).
 */
cmatrix_t *ucc_build_generator(int n_modes, const ucc_single_t *singles,
                               const double *theta_s, int n_singles,
                               const ucc_double_t *doubles,
                               const double *theta_d, int n_doubles);

/*
 * Prepare |\psi> = \exp(generator) |reference>: diagonalizes Hermitian matrix
 * i*generator (via cmatrix_eigh_complex) and applies \exp(generator) = V
 * \exp(-i * D) V^\dagger to `reference` without ever materializing dense
 * exp(generator) unitary itself. `generator` must be square and anti-Hermitian;
 * `reference` must match its dimension.
 *
 * Returns a newly-allocated, already-normalized cvector_t, or NULL on
 * invalid/mismatched input.
 */
cvector_t *ucc_prepare_state(const cmatrix_t *generator,
                             const cvector_t *reference);

/*
 * Build singles/doubles excitation lists (and \theta = converged CCSD amplitude
 * for each) directly from a ccsd_amplitudes_t (physics/ccsd.h, via
 * ccsd_run_ex).
 * Singles: every occ x virt pair (i,a) with theta_s = t1[i,a].
 * Doubles: every occ pair i<j and virt pair a<b (in amp->occ[]/amp->virt[]'s
 *          given ordering : i<j/a<b avoids double-counting same physical
 *          excitation via T2's antisymmetry, e.g. t2[i,j,a,b]=-t2[j,i,a,b])
 *          with \theta_d = t2[i,j,a,b].
 */
int ucc_excitations_from_ccsd_amplitudes(
    const ccsd_amplitudes_t *amp, ucc_single_t **singles_out,
    double **theta_s_out, int *n_singles_out, ucc_double_t **doubles_out,
    double **theta_d_out, int *n_doubles_out);

#endif
