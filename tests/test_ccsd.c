/*
 * Test: CCSD (Coupled Cluster Singles and Doubles), spin-orbital formulation,
 * built on molecular_hf.c's RHF orbitals and molecular_ao_to_mo's MO integrals.
 *
 * 1. H2/STO-3G: CCSD is exact for a 2-electron system (T1+T2 already span the
 *    full excitation space for only 2 electrons), so this checks against both
 *    an exact FCI diagonalization  (second_quant_build_molecular_hamiltonian)
 *    at the same geometry
 * 2. An asymmetric 4-hydrogen cluster (arbitrary non-collinear geometry, no
 *    spatial symmetry): a real >2-electron correlation test where CCSD is not
 *    exact.
 * 3. LiH/STO-3G, both with and without a frozen core: a real degenerate-orbital
 *    system (2px/2py exactly degenerate by symmetry).
 */

#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/ccsd.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/second_quant.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(const char *label, double got, double expected,
                        double tol) {
  double err = fabs(got - expected);
  printf("  %s: got=%.10f expected=%.10f err=%.2e\n", label, got, expected,
         err);
  if (err > tol) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAILED");
  if (!cond) {
    failures++;
  }
}

static void test_h2_ccsd_matches_fci(void) {
  printf("test_h2_ccsd_matches_fci:\n");
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, 1.4};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};
  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, 1.4}};
  molecule_t *mol = molecule_alloc(2, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 2, mol, 2, 1e-12, 200);

  check_true(res->converged, "H2 RHF converges");

  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 2, mol);
  double *eri_ao = molecular_eri_tensor(basis, 2);
  double *h_mo = malloc(4UL * sizeof(double));
  double *eri_mo = malloc(16UL * sizeof(double));

  molecular_ao_to_mo(h_ao, eri_ao, res->C, 2, h_mo, eri_mo);

  ccsd_result_t *ccsd = ccsd_run(2, h_mo, eri_mo, res->orbital_energies, 2, 0,
                                 res->total_energy, 1e-12, 100);
  check_true(ccsd->converged, "H2 CCSD converges");

  // independent reference 1: at same geometry
  check_close("H2 CCSD matches independent PySCF CCSD", ccsd->total_energy,
              -1.1372759436170439, 1e-8);

  // independent reference 2: exact FCI diagonalization - CCSD should be exact
  // here (2-electron system)
  cmatrix_t *H_fci =
      second_quant_build_molecular_hamiltonian(2, h_mo, eri_mo, 1.0 / 1.4);
  eigen_t *eig = cmatrix_eigh_complex(H_fci);
  double e_fci = eig->eigenvalues[0];

  for (int i = 1; i < eig->n; i++) {
    if (eig->eigenvalues[i] < e_fci) {
      e_fci = eig->eigenvalues[i];
    }
  }

  check_close("H2 CCSD matches this project's own exact FCI (exact for 2e-)",
              ccsd->total_energy, e_fci, 1e-8);

  eigen_free(eig);
  cmatrix_free(H_fci);
  free(ccsd);
  free(h_mo);
  free(eri_mo);
  free(eri_ao);
  cmatrix_free(h_ao);
  molecular_hf_result_free(res);
  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);
}

static void test_h4_asymmetric_ccsd_vs_pyscf(void) {
  printf("test_h4_asymmetric_ccsd_vs_pyscf:\n");
  // Arbitrary non-collinear, no-symmetry geometry
  double c0[3] = {0, 0, 0};
  double c1[3] = {0.9, 0.3, 0.1};
  double c2[3] = {1.7, -0.4, 0.6};
  double c3[3] = {2.9, 0.5, -0.3};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *h2 = molint_basis_sto3g_h(c2);
  basis_function_t *h3 = molint_basis_sto3g_h(c3);
  basis_function_t *basis[4] = {h0, h1, h2, h3};
  const double charges[4] = {1.0, 1.0, 1.0, 1.0};

  double centers[4][3];
  for (int d = 0; d < 3; d++) {
    centers[0][d] = c0[d];
    centers[1][d] = c1[d];
    centers[2][d] = c2[d];
    centers[3][d] = c3[d];
  }

  molecule_t *mol = molecule_alloc(4, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 4, mol, 4, 1e-12, 300);

  check_true(res->converged, "asymmetric H4 RHF converges");
  check_close("asymmetric H4 RHF", res->total_energy, -1.7475190952474993,
              1e-7);

  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 4, mol);
  double *eri_ao = molecular_eri_tensor(basis, 4);
  double *h_mo = malloc(16UL * sizeof(double));
  double *eri_mo = malloc(4UL * 4UL * 4UL * 4UL * sizeof(double));

  molecular_ao_to_mo(h_ao, eri_ao, res->C, 4, h_mo, eri_mo);

  ccsd_result_t *ccsd = ccsd_run(4, h_mo, eri_mo, res->orbital_energies, 4, 0,
                                 res->total_energy, 1e-12, 150);

  check_true(ccsd->converged, "asymmetric H4 CCSD converges");
  check_close("asymmetric H4 CCSD", ccsd->total_energy, -1.7894541829894628,
              1e-7);
  check_true(
      ccsd->total_energy < res->total_energy - 1e-4,
      "CCSD correlation energy is a real, non-trivial lowering below RHF");

  free(ccsd);
  free(h_mo);
  free(eri_mo);
  free(eri_ao);
  cmatrix_free(h_ao);
  molecular_hf_result_free(res);
  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);
  basis_function_free(h2);
  basis_function_free(h3);
}

static void test_lih_ccsd_vs_pyscf(void) {
  printf("test_lih_ccsd_vs_pyscf:\n");
  /* LiH/STO-3G, R=3.015 bohr: a real degenerate-orbital system (2px/2py
   * exactly degenerate by symmetry) */
  double cLi[3] = {0, 0, 0};
  double cH[3] = {0, 0, 3.015};
  basis_function_t *li[5];

  molint_basis_sto3g_li(cLi, li);

  basis_function_t *h = molint_basis_sto3g_h(cH);
  basis_function_t *basis[6] = {li[0], li[1], li[2], li[3], li[4], h};

  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, 3.015}};

  molecule_t *mol = molecule_alloc(2, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 6, mol, 4, 1e-12, 300);

  check_true(res->converged, "LiH RHF converges");
  check_close("LiH RHF", res->total_energy, -7.862009272120229, 1e-7);

  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 6, mol);
  double *eri_ao = molecular_eri_tensor(basis, 6);
  double *h_mo = malloc(36UL * sizeof(double));
  double *eri_mo = malloc(6UL * 6UL * 6UL * 6UL * sizeof(double));

  molecular_ao_to_mo(h_ao, eri_ao, res->C, 6, h_mo, eri_mo);

  ccsd_result_t *ccsd_full = ccsd_run(6, h_mo, eri_mo, res->orbital_energies, 4,
                                      0, res->total_energy, 1e-12, 200);

  check_true(ccsd_full->converged, "LiH (full) CCSD converges");
  check_close("LiH (full) CCSD matches independent PySCF CCSD",
              ccsd_full->total_energy, -7.882384455686598, 1e-6);

  ccsd_result_t *ccsd_frozen = ccsd_run(6, h_mo, eri_mo, res->orbital_energies,
                                        4, 1, res->total_energy, 1e-12, 200);

  check_true(ccsd_frozen->converged, "LiH (frozen-core) CCSD converges");
  check_close("LiH (frozen-core) CCSD frozen-core CCSD",
              ccsd_frozen->total_energy, -7.882167498191367, 1e-6);
  check_true(ccsd_frozen->total_energy > ccsd_full->total_energy - 1e-3,
             "frozen-core CCSD stays close to full-space CCSD (Li 1s core is "
             "nearly inert)");

  free(ccsd_full);
  free(ccsd_frozen);
  free(h_mo);
  free(eri_mo);
  free(eri_ao);
  cmatrix_free(h_ao);
  molecular_hf_result_free(res);
  molecule_free(mol);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }
  basis_function_free(h);
}

int main(void) {
  test_h2_ccsd_matches_fci();
  test_h4_asymmetric_ccsd_vs_pyscf();
  test_lih_ccsd_vs_pyscf();

  if (failures > 0) {
    printf("\n%d test_ccsd check(s) FAILED.\n", failures);
    return 1;
  }
  printf("\nAll test_ccsd checks passed.\n");
  return 0;
}
