/*
 * Test: CASSCF (Complete Active Space SCF)
 *
 * Reference values from mcscf.CASSCF() for LiH/STO-3G, R=3.015 bohr:
 *   CAS(2e,2o), 1 frozen core orbital: -7.881111895037824 Hartree
 *   CAS(4e,3o), no frozen core:        -7.881130996863274 Hartree
 */

#include "../physics/casscf.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
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
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

static void setup_lih(basis_function_t **basis, molecule_t **mol,
                      molecular_hf_result_t **hf) {
  double R = 3.015;
  static basis_function_t *li_orbs[5];
  double c_li[3] = {0, 0, 0}, c_h[3] = {0, 0, R};

  molint_basis_sto3g_li(c_li, li_orbs);
  basis_function_t *h_orb = molint_basis_sto3g_h(c_h);

  basis[0] = li_orbs[0];
  basis[1] = li_orbs[1];
  basis[2] = li_orbs[2];
  basis[3] = li_orbs[3];
  basis[4] = li_orbs[4];
  basis[5] = h_orb;

  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  *mol = molecule_alloc(2, charges, centers);

  *hf = molecular_rhf(basis, 6, *mol, 4, 1e-12, 200);
}

static void teardown_lih(basis_function_t **basis, molecule_t *mol,
                         molecular_hf_result_t *hf) {
  molecular_hf_result_free(hf);
  for (int i = 0; i < 5; i++) {
    basis_function_free(basis[i]);
  }
  basis_function_free(basis[5]);
  molecule_free(mol);
}

static void test_cas22_frozen_core(void) {
  printf("Test: LiH/STO-3G CASSCF(2e,2o), 1 frozen core orbital\n");

  basis_function_t *basis[6];
  molecule_t *mol;
  molecular_hf_result_t *hf;

  setup_lih(basis, &mol, &hf);

  casscf_result_t *r = casscf_run(basis, 6, mol, hf->C, 1, 2, 2, 1e-6, 60);

  check(r != NULL, "casscf_run should succeed");
  if (r) {
    check(r->converged, "CASSCF(2,2) converges within 60 iterations");
    check_close(r->total_energy, -7.881111895037824, 1e-5,
                "CASSCF(2,2) energy matches mcscf.CASSCF() reference");
    check(r->total_energy < hf->total_energy + 1e-6,
          "CASSCF energy is variationally at or below RHF");
    check(r->grad_norm < 1e-6,
          "final gradient norm below convergence threshold");

    casscf_result_free(r);
  }

  teardown_lih(basis, mol, hf);
}

static void test_cas43_no_frozen_core(void) {
  printf("Test: LiH/STO-3G CASSCF(4e,3o), no frozen core\n");

  basis_function_t *basis[6];
  molecule_t *mol;
  molecular_hf_result_t *hf;

  setup_lih(basis, &mol, &hf);

  casscf_result_t *r = casscf_run(basis, 6, mol, hf->C, 0, 3, 4, 1e-6, 60);

  check(r != NULL, "casscf_run should succeed");
  if (r) {
    check(r->converged, "CASSCF(4,3) converges within 60 iterations");
    check_close(r->total_energy, -7.881130996863274, 1e-5,
                "CASSCF(4,3) energy matches mcscf.CASSCF() reference");

    casscf_result_free(r);
  }

  teardown_lih(basis, mol, hf);
}

static void test_invalid_inputs_rejected(void) {
  printf("Test: invalid inputs are rejected cleanly\n");

  basis_function_t *basis[6];
  molecule_t *mol;
  molecular_hf_result_t *hf;

  setup_lih(basis, &mol, &hf);

  check(casscf_run(NULL, 6, mol, hf->C, 1, 2, 2, 1e-6, 60) == NULL,
        "NULL basis should be rejected");
  check(casscf_run(basis, 0, mol, hf->C, 1, 2, 2, 1e-6, 60) == NULL,
        "n_basis<1 should be rejected");
  check(casscf_run(basis, 6, mol, hf->C, -1, 2, 2, 1e-6, 60) == NULL,
        "negative n_frozen should be rejected");
  check(casscf_run(basis, 6, mol, hf->C, 1, 0, 2, 1e-6, 60) == NULL,
        "n_active<1 should be rejected");
  check(casscf_run(basis, 6, mol, hf->C, 5, 5, 2, 1e-6, 60) == NULL,
        "n_frozen+n_active>n_basis should be rejected");
  check(casscf_run(basis, 6, mol, hf->C, 1, 2, 5, 1e-6, 60) == NULL,
        "n_electrons_active > 2*n_active should be rejected");
  check(casscf_run(basis, 6, mol, hf->C, 1, 2, 2, 1e-6, 0) == NULL,
        "max_iter<1 should be rejected");

  teardown_lih(basis, mol, hf);
}

int main(void) {
  test_cas22_frozen_core();
  test_cas43_no_frozen_core();
  test_invalid_inputs_rejected();

  if (failures == 0) {
    printf("\nAll test_casscf checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
