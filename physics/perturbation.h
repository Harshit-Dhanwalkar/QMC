#ifndef QMC_PERTURBATION_H
#define QMC_PERTURBATION_H

#include "../core/matrix.h"
#include "../core/vector.h"

/* Non-degenerate perturbation theory (first and second order) */
typedef struct {
  double E0; // unperturbed energy
  double E1; // first-order correction
  double E2; // second-order correction
  /*
   * HACK: diagnostics for E2: perturb_nondeg's second-order sum omits any term
   * whose energy denominator |E_i - E_k| falls below `tol` (a near-degenerate
   * state k), since non-degenerate PT's 1/(E_i-E_k) formula diverges there.
   *
   * Where
   *   n_terms_skipped_near_degenerate: how many k were omitted (0 for a clean
   *                                    non-degenerate calculation).
   *   min_skipped_denominator        : min(|E_i-E_k|) among omitted terms (only
   *                                    meaningful when
   *                                    n_terms_skipped_near_degenerate > 0),
   *                                    closer to 0, more E2 likely
   * underestimates true correction.
   *
   * When n_terms_skipped_near_degenerate > 0, use perturb_degenerate on
   * near-degenerate subspace instead of trusting E2 from this call.
   */
  int n_terms_skipped_near_degenerate;
  double min_skipped_denominator;
} perturb_result_t;

perturb_result_t perturb_nondeg(const eigen_t *unperturbed, int state_index,
                                const cmatrix_t *V_pert, double tol);

/* Degenerate perturbation theory: diagonalize V_pert in the degenerate
 * subspace.
 * Returns new eigen_t with energies (eigenvalues) and eigenvectors (in that
 * subspace).
 * The subspace is defined by indices of states with same (or close) energy.
 *
 * Where
 *   energies          : array of unperturbed energies (size n)
 *   V_pert            : perturbation matrix (size n x n)
 *   degeneracy_indices: array of indices that belong to degenerate subspace.
 *
 * Returns eigen_t of the effective Hamiltonian (size = subspace size).
 */
eigen_t *perturb_degenerate(const double *energies, const cmatrix_t *V_pert,
                            const int *degeneracy_indices, int deg_size);

/* Time-dependent perturbation: Fermi's golden rule.
 * transition rate from initial state |i> to final states with density of states
 * \rho(E).
 *
 * Where
 *   V_pert  : perturbation operator (matrix element <f|V|i>)
 *   \rho_E  : density of states at energy E_f ~ E_i
 */
double fermi_golden_rate(const cmatrix_t *V_pert, int i, int f, double rho_E);

#endif
