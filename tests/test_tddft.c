/*
 * Test: Linear-response TDDFT, Tamm-Dancoff Approximation (TDA), for
 * closed-shell Kohn-Sham LDA (Slater exchange + PZ81 correlation)
 *
 * Reference values from tdscf.TDA(mf) with xc='LDA,PZ' (exactly this library's
 * LDA functional):
 *   H2/STO-3G:  1 excitation,  0.97123268 Hartree
 *   LiH/STO-3G: 5 excitations, 0.13833137, 0.18839063, 0.18839063 (doubly
 *               degenerate), 0.63802246, 1.73170728 Hartree
 * Both cross-checked to machine precision (1e-9 to 1e-14 Hartree) against
 * an independent implementation of same TDA formula using s own AO integrals
 * and numint xc-kernel machinery
 */

#include "../physics/molecular_dft.h"
#include "../physics/molecular_integrals.h"
#include "../physics/tddft.h"
#include <math.h>
#include <stdio.h>

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
    printf("  FAIL: %s (got %.8f, expected %.8f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

static void test_h2(void) {
  printf("  === Test: H2/STO-3G TDA-LDA (1 excitation) ===\n");

  double R = 1.4;
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};
  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);

  molecular_grid_t *grid = molecular_grid_build_default(mol);
  molecular_dft_result_t *ks = molecular_ks_lda_default(basis, 2, mol, 2, grid);

  check(ks != NULL && ks->converged, "H2 KS-LDA SCF converges");

  tddft_result_t *td = tddft_tda_lda(basis, 2, mol, ks, grid);
  check(td != NULL, "tddft_tda_lda should succeed");
  if (td) {
    check(td->dim == 1, "H2/STO-3G has exactly 1 possible excitation "
                        "(1 occ x 1 virt)");
    check_close(td->excitation_energies[0], 0.97123268, 1e-4,
                "H2/STO-3G TDA-LDA excitation energ matches referencey");

    tddft_result_free(td);
  }

  molecular_dft_result_free(ks);
  molecular_grid_free(grid);
  basis_function_free(h0);
  basis_function_free(h1);
  molecule_free(mol);
}

static void test_lih(void) {
  printf("  === Test: LiH/STO-3G TDA-LDA (5 excitations, incl. a degenerate "
         "pair) ===\n");

  double R = 3.015;
  double c_li[3] = {0, 0, 0}, c_h[3] = {0, 0, R};
  basis_function_t *li_orbs[5];
  molint_basis_sto3g_li(c_li, li_orbs);
  basis_function_t *h_orb = molint_basis_sto3g_h(c_h);
  basis_function_t *basis[6] = {li_orbs[0], li_orbs[1], li_orbs[2],
                                li_orbs[3], li_orbs[4], h_orb};
  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);

  molecular_grid_t *grid = molecular_grid_build_default(mol);
  molecular_dft_result_t *ks = molecular_ks_lda_default(basis, 6, mol, 4, grid);

  check(ks != NULL && ks->converged, "LiH KS-LDA SCF converges");

  tddft_result_t *td = tddft_tda_lda(basis, 6, mol, ks, grid);
  check(td != NULL, "tddft_tda_lda should succeed");
  if (td) {
    check(td->n_occ == 2 && td->n_virt == 4 && td->dim == 8,
          "LiH/STO-3G active space: 2 occ x 4 virt = 8-dim TDA problem");

    const double expected[5] = {0.13833137, 0.18839063, 0.18839063, 0.63802246,
                                1.73170728};
    for (int k = 0; k < 5; k++) {
      char label[96];

      snprintf(label, sizeof label,
               "LiH/STO-3G TDA-LDA excitation %d matches reference", k);
      check_close(td->excitation_energies[k], expected[k], 1e-4, label);
    }

    check(td->excitation_energies[0] <= td->excitation_energies[1] + 1e-12,
          "excitation energies are ascending");

    tddft_result_free(td);
  }

  molecular_dft_result_free(ks);
  molecular_grid_free(grid);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li_orbs[i]);
  }
  basis_function_free(h_orb);
  molecule_free(mol);
}

static void test_invalid_inputs(void) {
  printf("  === Test: invalid inputs are rejected cleanly ===\n");

  check(tddft_tda_lda(NULL, 2, NULL, NULL, NULL) == NULL,
        "all-NULL input should be rejected");
}

int main(void) {
  printf(" > TDDFT for H2 and LiH tests\n");

  test_h2();
  test_lih();
  test_invalid_inputs();

  if (failures == 0) {
    printf("\nAll test_tddft checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
