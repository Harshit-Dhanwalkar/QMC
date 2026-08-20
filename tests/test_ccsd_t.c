/*
 * Test: CCSD(T) perturbative triples correction, physics/ccsd_t.c.
 *
 * Reuses the exact same geometries as test_ccsd.c so RHF/CCSD reference values
 *
 * 1. H2/STO-3G (R=1.4 bohr): (T) is exactly zero here, not just small :
 *    forming a triples excitation needs 3 distinct occupied spin-orbital
 *    indices (i,j,k), but H2/STO-3G has only 2 occupied spin orbitals, so the
 *    triples excitation space is empty and every t3c/t3d term vanishes
 *    identically. This is a exact check on P(i/jk)P(a/bc) triple-occupied-index
 *    loop  bounds, not just a "close to zero" numerical check.
 * 2. Asymmetric H4 (same non-collinear, no-symmetry geometry as test_ccsd.c):
 *    (T) = -1.3256211e-05, a real nonzero correction with 4 occupied spin
 *    orbitals : exercises the full P(i/jk)P(a/bc) 9-term antisymmetrizer for
 *    both t3c and t3d.
 * 3. LiH/STO-3G (R=3.015 bohr, same as test_ccsd.c and test_lih.c): (T) =
 *    -8.393447e-06.
 */

#include "../core/matrix.h"
#include "../physics/ccsd_t.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
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

static void run_case(basis_function_t **basis, int n_basis, molecule_t *mol,
                     int n_electrons, double expected_rhf,
                     double expected_ccsd_corr, double expected_pert_t,
                     double tol, const char *label) {
  printf("test_ccsd_t_%s:\n", label);

  molecular_hf_result_t *hf =
      molecular_rhf(basis, n_basis, mol, n_electrons, 1e-12, 300);
  check_true(hf != NULL && hf->converged, "RHF converges");

  if (!hf) {
    return;
  }

  check_close(hf->total_energy, expected_rhf, 1e-7, "RHF matches reference");

  cmatrix_t *Hcore = molecular_core_hamiltonian(basis, n_basis, mol);
  double *eri_ao = molecular_eri_tensor(basis, n_basis);
  double *h_mo = malloc((size_t)n_basis * n_basis * sizeof(double));
  double *eri_mo =
      malloc((size_t)n_basis * n_basis * n_basis * n_basis * sizeof(double));

  molecular_ao_to_mo(Hcore, eri_ao, hf->C, n_basis, h_mo, eri_mo);

  ccsdt_result_t *ccsdt =
      ccsdt_run(n_basis, h_mo, eri_mo, hf->orbital_energies, n_electrons, 0,
                hf->total_energy, 1e-12, 200);
  check_true(ccsdt != NULL, "ccsdt_run allocates");

  if (ccsdt) {
    check_true(ccsdt->ccsd_converged, "CCSD converges");
    printf("  CCSD corr=%.10f  (T)=%.10f  total=%.10f\n",
           ccsdt->ccsd_correlation_energy, ccsdt->perturbative_correction,
           ccsdt->total_energy);

    check_close(ccsdt->ccsd_correlation_energy, expected_ccsd_corr, tol,
                "CCSD correlation energy matches PySCF reference");
    check_close(ccsdt->perturbative_correction, expected_pert_t, tol,
                "(T) correction matches PySCF reference");
    check_close(ccsdt->total_energy,
                expected_rhf + expected_ccsd_corr + expected_pert_t, tol,
                "CCSD(T) total energy self-consistent");

    free(ccsdt);
  }

  free(h_mo);
  free(eri_mo);
  free(eri_ao);
  cmatrix_free(Hcore);
  molecular_hf_result_free(hf);
}

static void test_h2_pert_t_exactly_zero(void) {
  double R = 1.4;
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};
  const double charge[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charge, centers);

  run_case(basis, 2, mol, 2, -1.116714325063, -0.020561618563, 0.0, 1e-9,
           "h2_zero");

  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);
}

static void test_asymmetric_h4(void) {
  double c0[3] = {0, 0, 0};
  double c1[3] = {0.9, 0.3, 0.1};
  double c2[3] = {1.7, -0.4, 0.6};
  double c3[3] = {2.9, 0.5, -0.3};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *h2 = molint_basis_sto3g_h(c2);
  basis_function_t *h3 = molint_basis_sto3g_h(c3);
  basis_function_t *basis[4] = {h0, h1, h2, h3};
  const double charge[4] = {1.0, 1.0, 1.0, 1.0};
  double centers[4][3];
  for (int d = 0; d < 3; d++) {
    centers[0][d] = c0[d];
    centers[1][d] = c1[d];
    centers[2][d] = c2[d];
    centers[3][d] = c3[d];
  }

  molecule_t *mol = molecule_alloc(4, charge, centers);

  run_case(basis, 4, mol, 4, -1.747519095248, -0.041935087838, -0.000013256211,
           1e-7, "asymmetric_h4");

  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);
  basis_function_free(h2);
  basis_function_free(h3);
}

static void test_lih(void) {
  double R = 3.015;
  double li_center[3] = {0, 0, 0};
  basis_function_t *li[5];
  molint_basis_sto3g_li(li_center, li);
  double h_center[3] = {0, 0, R};
  basis_function_t *h = molint_basis_sto3g_h(h_center);
  basis_function_t *basis[6] = {li[0], li[1], li[2], li[3], li[4], h};
  const double charge[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charge, centers);

  run_case(basis, 6, mol, 4, -7.862009272120, -0.020375183567, -0.000008393447,
           1e-7, "lih");

  molecule_free(mol);
  basis_function_free(h);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }
}

int main(void) {
  test_h2_pert_t_exactly_zero();
  test_asymmetric_h4();
  test_lih();

  if (failures == 0) {
    printf("\nAll CCSD(T) tests PASSED\n");
    return 0;
  }
  printf("\n%d CCSD(T) test(s) FAILED\n", failures);

  return 1;
}
