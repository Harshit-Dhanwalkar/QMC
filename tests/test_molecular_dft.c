/*
 * Test: molecular (multi-atom) Kohn-Sham LDA DFT via Becke fuzzy-cell grids.
 *
 *   1. Overlap-matrix reproduction: integrate S numerically on Becke grid
 *      (\sum_g weight_g * \phi_p(r_g) * \phi_q(r_g)) and compare against exact
 *      analytic overlap matrix from molecular_integrals.c. This isolates the
 *      grid construction (Becke weights, radial/angular quadrature) from SCF
 *      machinery entirely.
 *   2. Electron-count integration: sum_g weight_g * n(r_g) at a known density
 *      matrix must equal the exact electron count.
 *   3. End-to-end total energy on H2/STO-3G and LiH/STO-3G, cross-validated
 *      (Reference: Perdew-Zunger 1981 correlation), lda_xc_energy_density /
 *      lda_xc_potential in dft.c): H2/STO-3G @ R=1.4 bohr: -1.12132825509958
 *      Hartree LiH/STO-3G @ R=3.015 bohr: -7.79120636378942 Hartree
 */

#include "../physics/molecular_dft.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../core/matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

static void test_grid_reproduces_overlap_matrix(void) {
  printf("Test: Becke grid numerically reproduces the exact analytic "
         "overlap matrix (H2/STO-3G)\n");

  double R = 1.4;
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};
  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);
  cmatrix_t *S_exact = molecular_overlap_matrix(basis, 2);

  molecular_grid_t *grid = molecular_grid_build(mol, 60, 20, 30, 1.0);
  check(grid != NULL, "grid should build successfully");

  if (grid) {
    double S_num[2][2] = {{0, 0}, {0, 0}};

    for (int g = 0; g < grid->n_points; g++) {
      double r[3] = {grid->points[g].x, grid->points[g].y, grid->points[g].z};
      double phi0 = basis_function_value(h0, r);
      double phi1 = basis_function_value(h1, r);
      double w = grid->points[g].weight;

      S_num[0][0] += w * phi0 * phi0;
      S_num[0][1] += w * phi0 * phi1;
      S_num[1][1] += w * phi1 * phi1;
    }

    check_close(S_num[0][0], CMAT(S_exact, 0, 0).re, 1e-5,
                "numeric S[0][0] matches analytic overlap");
    check_close(S_num[0][1], CMAT(S_exact, 0, 1).re, 1e-5,
                "numeric S[0][1] matches analytic overlap");
    check_close(S_num[1][1], CMAT(S_exact, 1, 1).re, 1e-5,
                "numeric S[1][1] matches analytic overlap");
  }

  molecular_grid_free(grid);
  cmatrix_free(S_exact);
  basis_function_free(h0);
  basis_function_free(h1);
  molecule_free(mol);
}

static void test_grid_integrates_electron_count(void) {
  printf("Test: Becke grid correctly integrates a known density to the exact "
         "electron count (LiH/STO-3G, RHF density matrix)\n");

  double R = 3.015;
  basis_function_t *li_orbs[5];
  double c_li[3] = {0, 0, 0};
  double c_h[3] = {0, 0, R};
  molint_basis_sto3g_li(c_li, li_orbs);

  basis_function_t *h_orb = molint_basis_sto3g_h(c_h);
  basis_function_t *basis[6] = {li_orbs[0], li_orbs[1], li_orbs[2],
                                li_orbs[3], li_orbs[4], h_orb};
  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);

  molecular_hf_result_t *hf = molecular_rhf(basis, 6, mol, 4, 1e-10, 200);
  check(hf != NULL && hf->converged,
        "RHF should converge for the density matrix construction");

  molecular_grid_t *grid = molecular_grid_build(mol, 60, 20, 30, 1.0);
  check(grid != NULL, "grid should build successfully");

  if (hf && grid) {
    // Build the RHF density matrix from occupied MOs
    double D[6][6];

    for (int p = 0; p < 6; p++) {
      for (int q = 0; q < 6; q++) {
        double v = 0.0;

        for (int k = 0; k < 2; k++) { /* n_occ = 4 electrons / 2 = 2 */
          v += 2.0 * CMAT(hf->C, p, k).re * CMAT(hf->C, q, k).re;
        }

        D[p][q] = v;
      }
    }

    double n_electrons_numeric = 0.0;
    for (int g = 0; g < grid->n_points; g++) {
      double r[3] = {grid->points[g].x, grid->points[g].y, grid->points[g].z};
      double phi[6];

      for (int p = 0; p < 6; p++) {
        phi[p] = basis_function_value(basis[p], r);
      }

      double dens = 0.0;

      for (int p = 0; p < 6; p++) {
        for (int q = 0; q < 6; q++) {
          dens += D[p][q] * phi[p] * phi[q];
        }
      }

      n_electrons_numeric += grid->points[g].weight * dens;
    }

    check_close(n_electrons_numeric, 4.0, 1e-4,
                "grid-integrated electron count matches exact N=4");
  }

  molecular_grid_free(grid);
  if (hf) {
    molecular_hf_result_free(hf);
  }
  for (int i = 0; i < 5; i++) {
    basis_function_free(li_orbs[i]);
  }
  basis_function_free(h_orb);
  molecule_free(mol);
}

static void test_h2_ks_lda_matches_pyscf(void) {
  printf("Test: H2/STO-3G KS-LDA total energy matches PySCF's "
         "dft.RKS(xc='lda,pz') reference (-1.12132825509958 Hartree)\n");

  double R = 1.4;
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};
  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);

  molecular_grid_t *grid = molecular_grid_build_default(mol);
  check(grid != NULL, "grid should build");

  if (grid) {
    molecular_dft_result_t *dft =
        molecular_ks_lda_default(basis, 2, mol, 2, grid);
    check(dft != NULL, "KS-LDA should allocate a result");

    if (dft) {
      check(dft->converged, "KS-LDA should converge");
      check_close(dft->total_energy, -1.12132825509958, 1e-5,
                  "H2/STO-3G KS-LDA total energy matches PySCF reference");

      molecular_dft_result_free(dft);
    }
  }

  molecular_grid_free(grid);
  basis_function_free(h0);
  basis_function_free(h1);
  molecule_free(mol);
}

static void test_lih_ks_lda_matches_pyscf(void) {
  printf("Test: LiH/STO-3G KS-LDA total energy matches PySCF's "
         "dft.RKS(xc='lda,pz') reference (-7.79120636378942 Hartree), "
         "requires density-mixing convergence\n");

  double R = 3.015;
  basis_function_t *li_orbs[5];
  double c_li[3] = {0, 0, 0};
  double c_h[3] = {0, 0, R};

  molint_basis_sto3g_li(c_li, li_orbs);
  basis_function_t *h_orb = molint_basis_sto3g_h(c_h);
  basis_function_t *basis[6] = {li_orbs[0], li_orbs[1], li_orbs[2],
                                li_orbs[3], li_orbs[4], h_orb};

  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);

  molecular_grid_t *grid = molecular_grid_build_default(mol);
  check(grid != NULL, "grid should build");

  if (grid) {
    molecular_dft_result_t *dft =
        molecular_ks_lda_default(basis, 6, mol, 4, grid);
    check(dft != NULL, "KS-LDA should allocate a result");

    if (dft) {
      check(dft->converged, "KS-LDA should converge (with mixing)");
      check_close(dft->total_energy, -7.79120636378942, 1e-5,
                  "LiH/STO-3G KS-LDA total energy matches PySCF reference");

      check_close(dft->e_core + dft->e_coulomb + dft->e_xc + dft->e_nuclear,
                  dft->total_energy, 1e-4,
                  "reported energy components sum to the total energy");

      molecular_dft_result_free(dft);
    }
  }

  molecular_grid_free(grid);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li_orbs[i]);
  }
  basis_function_free(h_orb);
  molecule_free(mol);
}

static void test_invalid_inputs_rejected(void) {
  printf("Test: invalid inputs are rejected cleanly (NULL, not crashes)\n");

  check(molecular_grid_build(NULL, 10, 10, 10, 1.0) == NULL,
        "NULL molecule should be rejected");

  double c0[3] = {0, 0, 0};
  const double charges[1] = {1.0};
  double centers[1][3] = {{0, 0, 0}};
  molecule_t *mol = molecule_alloc(1, charges, centers);

  check(molecular_grid_build(mol, 0, 10, 10, 1.0) == NULL,
        "zero radial points should be rejected");
  check(molecular_grid_build(mol, 10, 10, 10, -1.0) == NULL,
        "negative radial scale should be rejected");

  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *basis[1] = {h0};
  molecular_grid_t *grid = molecular_grid_build_default(mol);

  check(molecular_ks_lda(basis, 1, mol, 1, grid, 0.3, 1e-9, 100) == NULL,
        "odd electron count should be rejected (restricted closed-shell only)");
  check(molecular_ks_lda(NULL, 1, mol, 2, grid, 0.3, 1e-9, 100) == NULL,
        "NULL basis should be rejected");

  molecular_grid_free(grid);
  basis_function_free(h0);
  molecule_free(mol);
}

int main(void) {
  printf("=== Molecular Kohn-Sham LDA DFT (Becke grid) tests ===\n\n");

  test_grid_reproduces_overlap_matrix();
  test_grid_integrates_electron_count();
  test_h2_ks_lda_matches_pyscf();
  test_lih_ks_lda_matches_pyscf();
  test_invalid_inputs_rejected();

  if (failures == 0) {
    printf("\nAll test_molecular_dft checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_molecular_dhf check(s) FAILED.\n", failures);
    return 1;
  }
}
