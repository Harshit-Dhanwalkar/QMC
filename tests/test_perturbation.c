/*
 * Test perturbation theory: compare first-order correction for harmonic + x^4
 * with known analytic result for the ground state.
 * For atomic units m = \omega=1, <0|x^4|0> = 3/4.
 *
 * Validates
 *  1. perturb_nondeg's first- and second-order corrections against
 *     harmonic-oscillator + \lambda * x^4 system,
 *  2. perturb_nondeg's near-degeneracy diagnostic
 *     (n_terms_skipped_near_degenerate, min_skipped_denominator) on a toy
 *     2-level system engineered to trigger it, cross-validated against
 *     perturb_degenerate's exact splitting for that same pair,
 *  3. fermi_golden_rate against its closed-form definition directly
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

#ifndef RUNNING_ON_VALGRIND
#define RUNNING_ON_VALGRIND 0
#endif

static int failures = 0;

static void check(int cond, const char *msg) {
  if (!cond) {
    printf("  FAIL: %s\n", msg);
    failures++;
  }
}

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  if (fabs(got - expected) > tol) {
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

static double eigenbasis_element(const eigen_t *eig, int row, int col,
                                 const double *V_grid, int N, double dx) {
  double norm_row = 0.0, norm_col = 0.0;
  double element = 0.0;
  for (int i = 0; i < N; i++) {
    double a = CMAT(eig->eigenvectors, i, row).re;
    double b = CMAT(eig->eigenvectors, i, col).re;

    norm_row += a * a;
    norm_col += b * b;
    element += a * V_grid[i] * b;
  }

  norm_row = sqrt(norm_row * dx);
  norm_col = sqrt(norm_col * dx);

  return element * dx / (norm_row * norm_col);
}

/* Test 1: perturb_nondeg on harmonic oscillator + \lambda*x^4, cross-checked
 * against full numerical diagonalization */
static void test_perturb_nondeg_anharmonic_oscillator(void) {
  printf("Test: perturb_nondeg on harmonic + lambda*x^4 vs. full "
         "diagonalization\n");

  double omega = 1.0, m = 1.0, hbar = 1.0, lambda = 0.05;
  int N = RUNNING_ON_VALGRIND ? 101 : 501;
  double x_min = -6.0, x_max = 6.0;
  double dx = (x_max - x_min) / (N - 1);
  double *x = linspace(x_min, x_max, N);

  check(x != NULL, "linspace should succeed");
  if (!x) {
    return;
  }

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
  check(eig0 != NULL, "H_0 diagonalization should succeed");

  cmatrix_t *H_full = cmatrix_copy(H_0);
  for (int i = 0; i < N; i++) {
    CMAT(H_full, i, i) = c_add(CMAT(H_full, i, i), c_real(V_grid[i]));
  }

  eigen_t *eig_full = cmatrix_eigh(H_full);
  check(eig_full != NULL, "full H diagonalization should succeed");

  if (eig0 && eig_full) {
    // Build only the eigenbasis columns needed (states 0..2)
    int n_states = 3;
    cmatrix_t *V_pert = cmatrix_alloc(N, N);

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        CMAT(V_pert, i, j) = c_zero();
      }
    }

    for (int col = 0; col < n_states; col++) {
      for (int k = 0; k < N; k++) {
        double elt = eigenbasis_element(eig0, k, col, V_grid, N, dx);
        CMAT(V_pert, k, col) = c_real(elt);
        CMAT(V_pert, col, k) = c_real(elt);
      }
    }

    /* n=0 first-order correction has an exact closed form:
     * E1 = \lambda * <0|x^4|0> = \lambda * (3/4) for m = \hbar = \omega = 1 QHO
     */
    perturb_result_t res0 = perturb_nondeg(eig0, 0, V_pert, 1e-8);
    check_close(res0.E1, lambda * 0.75, 1e-3,
                "n=0 first-order correction matches exact <0|x^4|0>=3/4");
    check(res0.n_terms_skipped_near_degenerate == 0,
          "no near-degenerate terms should be skipped for a healthy "
          "non-degenerate oscillator spectrum");
    /* Full E0 + E1 + E2 vs numerical diagonalization, for n=0,1,2. Second-order
     * PT for a weak (\lambda=0.05) perturbation should land within a few
     * percent of the numerically exact value at these low levels. */
    for (int n = 0; n < n_states; n++) {
      perturb_result_t res = perturb_nondeg(eig0, n, V_pert, 1e-8);

      double E_pt = res.E0 + res.E1 + res.E2;
      double E_full = eig_full->eigenvalues[n];
      double rel_err = fabs(E_full - E_pt) / fabs(E_full);

      char msg[128];
      snprintf(
          msg, sizeof(msg),
          "n=%d: perturb_nondeg's E0+E1+E2 within 3%% of full diagonalization",
          n);

      check(rel_err < 0.03, msg);
    }

    cmatrix_free(V_pert);
  }

  cmatrix_free(H_0);
  cmatrix_free(H_full);
  if (eig0)
    eigen_free(eig0);
  if (eig_full)
    eigen_free(eig_full);
  free(V_grid);
  free(x);
}

/* Test 2: near-degeneracy diagnostic, cross-validated against
 * perturb_degenerate's exact splitting on the same pair. */
static void test_perturb_nondeg_near_degeneracy_diagnostic(void) {
  printf("Test: perturb_nondeg flags near-degenerate perturb_degenerate "
         "gives splitting for same pair\n");

  double gap = 1e-9; // far below any reasonable tol
  double coupling = 0.01;
  double energies[2] = {1.0, 1.0 + gap};
  eigen_t eig = {2, energies, NULL};

  cmatrix_t *V = cmatrix_alloc(2, 2);
  CMAT(V, 0, 0) = c_real(0.0);
  CMAT(V, 1, 1) = c_real(0.0);
  CMAT(V, 0, 1) = c_real(coupling);
  CMAT(V, 1, 0) = c_real(coupling);

  double tol = 1e-6; // >> gap, so the k=1 term must be flagged for state 0
  perturb_result_t res = perturb_nondeg(&eig, 0, V, tol);

  check(res.n_terms_skipped_near_degenerate == 1,
        "exactly one near-degenerate term should be flagged");
  check_close(res.min_skipped_denominator, gap, 1e-12,
              "reported skipped denominator matches true (tiny) energy gap");
  check_close(res.E2, 0.0, 1e-12, "E2 should not include near-degenerate term");

  /* Cross-check: perturb_degenerate on the same 2-level subspace gives exact,
   * splitting E +/- |V_01| (to O(gap), negligible here) */
  const int deg_indices[2] = {0, 1};
  eigen_t *deg_eig = perturb_degenerate(energies, V, deg_indices, 2);
  check(deg_eig != NULL, "perturb_degenerate should succeed");
  if (deg_eig) {
    double lo = fmin(deg_eig->eigenvalues[0], deg_eig->eigenvalues[1]);
    double hi = fmax(deg_eig->eigenvalues[0], deg_eig->eigenvalues[1]);
    check_close(hi - lo, 2.0 * coupling, 1e-6,
                "perturb_degenerate recovers exact splitting 2*|V_01| that "
                "perturb_nondeg correctly declines to compute");

    eigen_free(deg_eig);
  }

  cmatrix_free(V);
}

/* Test 3: fermi_golden_rate against its closed-form definition
 *  \Gamma = 2 * \pi * |V_fi|^2 * \rho(E), directly. */
static void test_fermi_golden_rate(void) {
  printf("Test: fermi_golden_rate matches Gamma=2*pi*|V_fi|^2*rho(E) "
         "directly\n");

  cmatrix_t *V = cmatrix_alloc(2, 2);
  CMAT(V, 0, 1) = c_new(0.3, 0.4); // |V_fi|^2 = 0.25
  CMAT(V, 1, 0) = c_new(0.3, -0.4);

  double rho_E = 2.0;
  double rate = fermi_golden_rate(V, 0, 1, rho_E);
  double expected = 2.0 * M_PI * 0.25 * rho_E;

  check_close(rate, expected, 1e-12,
              "fermi_golden_rate matches the closed-form Fermi's Golden "
              "Rule formula exactly");
  check_close(fermi_golden_rate(NULL, 0, 1, rho_E), 0.0, 1e-15,
              "NULL V_pert returns 0 rather than crashing");
  check_close(fermi_golden_rate(V, -1, 1, rho_E), 0.0, 1e-15,
              "negative index returns 0 rather than an out-of-bounds read");

  cmatrix_free(V);
}

static void test_ground_state_x4_expectation(void) {
  printf("Test: ground-state <x^4> matches exact analytic result 3/4\n");

  int N = RUNNING_ON_VALGRIND ? 51 : 501;
  double x_min = -6.0, x_max = 6.0;
  double dx = (x_max - x_min) / (N - 1);
  double *x = linspace(x_min, x_max, N);

  check(x != NULL, "linspace should succeed");
  if (!x) {
    return;
  }

  double coeff = -0.5 / (dx * dx);
  cmatrix_t *H_0 = cmatrix_alloc(N, N);
  for (int i = 0; i < N; i++) {
    double V0 = 0.5 * x[i] * x[i];
    CMAT(H_0, i, i) = c_real(-2.0 * coeff + V0);

    if (i > 0) {
      CMAT(H_0, i, i - 1) = c_real(coeff);
    }
    if (i < N - 1) {
      CMAT(H_0, i, i + 1) = c_real(coeff);
    }
  }

  eigen_t *eig = cmatrix_eigh(H_0);
  check(eig != NULL, "diagonalization should succeed");
  if (eig) {
    check_close(eig->eigenvalues[0], 0.5, 0.05,
                "ground state energy close to exact 0.5");

    double norm_sq = 0.0;
    for (int i = 0; i < N; i++) {
      norm_sq += pow(CMAT(eig->eigenvectors, i, 0).re, 2);
    }
    norm_sq *= dx;

    double x4_expect = 0.0;
    for (int i = 0; i < N; i++) {
      double psi_i = CMAT(eig->eigenvectors, i, 0).re / sqrt(norm_sq);
      double xi2 = x[i] * x[i];
      x4_expect += xi2 * xi2 * psi_i * psi_i;
    }
    x4_expect *= dx;

    check_close(x4_expect, 0.75, 0.02,
                "<0|x^4|0> matches exact analytic value 3/4");

    eigen_free(eig);
  }

  cmatrix_free(H_0);
  free(x);
}

int main(void) {
  printf("=== physics/perturbation.h tests ===\n\n");

  test_perturb_nondeg_anharmonic_oscillator();
  test_perturb_nondeg_near_degeneracy_diagnostic();
  test_fermi_golden_rate();
  test_ground_state_x4_expectation();

  if (failures == 0) {
    printf("\nAll tests passed.\n");
    return 0;
  } else {
    printf("\n%d test(s) FAILED.\n", failures);
    return 1;
  }
}
