/*
1st/2nd order, degenerate, Stark/Zeeman
Perturbation theory:
Time‑independent (non‑deg/deg) and time‑dependent (Fermi golden rule)
*/

#include "perturbation.h"
#include "../core/complex.h"
#include "../core/linalg/linalg.h"
#include "../core/matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

perturb_result_t perturb_nondeg(const eigen_t *unperturbed, int state_index,
                                const cmatrix_t *V_pert, double tol) {
  perturb_result_t res = {0.0, 0.0, 0.0, 0, 0.0};
  if (!unperturbed || !V_pert || state_index < 0 ||
      state_index >= unperturbed->n) {
    return res;
  }

  int n = unperturbed->n;
  double E0 = unperturbed->eigenvalues[state_index];
  res.E0 = E0;

  // First-order: <i|V|i>
  complex_t Vii = CMAT(V_pert, state_index, state_index);
  res.E1 = Vii.re; // assume real Hermitian

  // Second-order: \sum_{k != i} |<k|V|i>|^2 / (E_i - E_k)
  double E2 = 0.0;
  int n_skipped = 0;
  double min_skipped_denom = 0.0;
  for (int k = 0; k < n; k++) {
    if (k == state_index) {
      continue;
    }

    double denom = E0 - unperturbed->eigenvalues[k];
    if (fabs(denom) < tol) {
      // silently under-counting.
      if (n_skipped == 0 || fabs(denom) < min_skipped_denom) {
        min_skipped_denom = fabs(denom);
      }

      n_skipped++;
      continue;
    }

    complex_t Vki = CMAT(V_pert, k, state_index);
    double num = c_abs2(Vki);

    E2 += num / denom;
  }

  res.E2 = E2;
  res.n_terms_skipped_near_degenerate = n_skipped;
  res.min_skipped_denominator = min_skipped_denom;

  return res;
}

eigen_t *perturb_degenerate(const double *energies, const cmatrix_t *V_pert,
                            const int *deg_indices, int deg_size) {
  if (!energies || !V_pert || !deg_indices || deg_size < 1) {
    return NULL;
  }

  // Build effective Hamiltonian (subspace)
  cmatrix_t *H_eff = cmatrix_alloc(deg_size, deg_size);
  if (!H_eff) {
    return NULL;
  }

  for (int i = 0; i < deg_size; i++) {
    int ii = deg_indices[i];

    for (int j = 0; j < deg_size; j++) {
      int jj = deg_indices[j];
      CMAT(H_eff, i, j) = CMAT(V_pert, ii, jj);
    }
  }

  eigen_t *eig = cmatrix_eigh_generic(H_eff);
  cmatrix_free(H_eff);

  return eig;
}

double fermi_golden_rate(const cmatrix_t *V_pert, int i, int f, double rho_E) {
  if (!V_pert || i < 0 || f < 0) {
    return 0.0;
  }

  complex_t Vfi = CMAT(V_pert, f, i);
  double rate = 2.0 * M_PI * c_abs2(Vfi) * rho_E;

  return rate;
}
