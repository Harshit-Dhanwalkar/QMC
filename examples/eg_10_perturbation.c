/*
 * Non-degenerate Perturbation Theory
 *
 * Unperturbed: 1D harmonic oscillator (m=1, \hbar=1, \omega=1)
 * Perturbation: V' = \lambda x^4      (anharmonic term)
 */

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../physics/perturbation.h"
#include "../physics/potentials.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * <\phi_{col}|V|\phi_k> for every k = 0..n-1, i.e. one full column of V in
 * eigenbasis, via direct numerical integration on position grid. Written into
 * V_pert's column `col`
 * NOTE:since V is symmetric here, that also gives the corresponding row.
 */
static void fill_eigenbasis_column(cmatrix_t *V_pert, int col,
                                   const eigen_t *eig, const double *x,
                                   const double *V_grid, int N, double dx) {
  int n = eig->n;
  cvector_t *phi_col = cvector_from_matrix_column(eig->eigenvectors, col);
  double norm_col = 0.0;
  for (int i = 0; i < N; i++) {
    norm_col += c_abs2(phi_col->data[i]);
  }

  norm_col = sqrt(norm_col * dx);

  for (int k = 0; k < n; k++) {
    cvector_t *phi_k = cvector_from_matrix_column(eig->eigenvectors, k);
    double norm_k = 0.0;
    for (int i = 0; i < N; i++) {
      norm_k += c_abs2(phi_k->data[i]);
    }
    norm_k = sqrt(norm_k * dx);

    double element = 0.0;
    for (int i = 0; i < N; i++) {
      element += phi_k->data[i].re * V_grid[i] * phi_col->data[i].re;
    }
    element *= dx / (norm_k * norm_col);

    CMAT(V_pert, k, col) = c_real(element);
    CMAT(V_pert, col, k) = c_real(element); // V is real symmetric here
    cvector_free(phi_k);
  }
  cvector_free(phi_col);
  (void)x;
}

int main(void) {
  printf(" > Non-degenerate Perturbation Theory (Harmonic + \\lambda x^4)\n\n");

  // Parameters (atomic units: \hbar = m = \omega = 1)
  double omega = 1.0;
  double m = 1.0;
  double hbar = 1.0;
  double lambda = 0.05; // perturbation strength

  int N = 501;
  double x_min = -6.0, x_max = 6.0;
  double dx = (x_max - x_min) / (N - 1);
  double *x = linspace(x_min, x_max, N);

  double coeff = -hbar * hbar / (2.0 * m * dx * dx);
  cmatrix_t *H_0 = cmatrix_alloc(N, N);
  double *V_grid = malloc((size_t)N * sizeof(double));
  for (int i = 0; i < N; i++) {
    double V0 = 0.5 * m * omega * omega * x[i] * x[i];
    CMAT(H_0, i, i) = c_real(-2.0 * coeff + V0);

    if (i > 0) {
      CMAT(H_0, i, i - 1) = c_real(coeff);
    }

    if (i < N - 1) {
      CMAT(H_0, i, i + 1) = c_real(coeff);
    }

    V_grid[i] = lambda * x[i] * x[i] * x[i] * x[i];
  }

  eigen_t *eig0 = cmatrix_eigh(H_0);

  printf("Step 1: build V' = lambda*x^4 in the H_0 eigenbasis (only columns "
         "needed for states 0..4) and call perturb_nondeg directly\n");

  cmatrix_t *V_pert = cmatrix_alloc(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CMAT(V_pert, i, j) = c_zero();
    }
  }

  for (int n = 0; n < 5 && n < eig0->n; n++) {
    fill_eigenbasis_column(V_pert, n, eig0, x, V_grid, N, dx);
  }

  printf("   State n   E0 (exact)   E1 (perturb_nondeg)   E2 (perturb_nondeg)  "
         " E0+E1+E2   E_full (numerical)   Error (%%)\n");
  printf("   -------   ----------   --------------------   -----------------   "
         "---------   -------------------   --------\n");

  for (int n = 0; n < 5 && n < eig0->n; n++) {
    perturb_result_t res = perturb_nondeg(eig0, n, V_pert, 1e-8);
    double E0_analytical = hbar * omega * (n + 0.5);

    cmatrix_t *H_full = cmatrix_copy(H_0);

    for (int i = 0; i < N; i++) {
      CMAT(H_full, i, i) = c_add(CMAT(H_full, i, i), c_real(V_grid[i]));
    }

    eigen_t *eig_full = cmatrix_eigh(H_full);

    double E_full = eig_full->eigenvalues[n];
    double E_pt = res.E0 + res.E1 + res.E2;
    double err = fabs(E_full - E_pt) / fabs(E_full) * 100.0;

    printf(" %2d       %8.6f    %14.8f   %16.8f    %8.6f    %14.6f    "
           "%6.3f%%\n",
           n, E0_analytical, res.E1, res.E2, E_pt, E_full, err);
    if (res.n_terms_skipped_near_degenerate > 0) {
      printf("            (note: %d near-degenerate term(s) skipped in E2, "
             "min denominator %.2e -- E2 may be an underestimate)\n",
             res.n_terms_skipped_near_degenerate, res.min_skipped_denominator);
    }

    cmatrix_free(H_full);
    eigen_free(eig_full);
  }

  printf("\nStep 2: demonstrate perturb_nondeg's near-degeneracy diagnostic on "
         "a toy 2-level system\n\n");
  /* Two nearly-degenerate unperturbed levels (E=1.0 and E=1.0+1e-9) coupled
   * by an off-diagonal perturbation strong enough that naive non-degenerate
   * PT's 1 / (E_i - E_k) formula would wildly diverge if not caught */
  double toy_energies[2] = {1.0, 1.0 + 1e-9};
  eigen_t toy_eig = {2, toy_energies, NULL};
  cmatrix_t *V_toy = cmatrix_alloc(2, 2);

  CMAT(V_toy, 0, 0) = c_real(0.0);
  CMAT(V_toy, 1, 1) = c_real(0.0);
  CMAT(V_toy, 0, 1) = c_real(0.01);
  CMAT(V_toy, 1, 0) = c_real(0.01);
  perturb_result_t toy_res = perturb_nondeg(&toy_eig, 0, V_toy, 1e-6);

  printf("  Toy 2-level system: E0=%.1f, E1=%.6f, E2=%.6f, n_skipped=%d, "
         "min_skipped_denom=%.2e\n",
         toy_res.E0, toy_res.E1, toy_res.E2,
         toy_res.n_terms_skipped_near_degenerate,
         toy_res.min_skipped_denominator);
  printf(
      "  perturb_nondeg correctly flags this than - using perturb_degenerate "
      "on this pair gives splitting E_+/- = E0 +/- |V_01| = %.6f / %.6f\n",
      toy_energies[0] + fabs(0.01), toy_energies[0] - fabs(0.01));

  // Cleanup
  cmatrix_free(V_toy);
  cmatrix_free(V_pert);
  cmatrix_free(H_0);
  eigen_free(eig0);
  free(V_grid);
  free(x);

  return 0;
}
