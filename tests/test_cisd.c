/*
 * Test: Configuration Interaction Singles and Doubles (CISD)
 *
 * Validation:
 *   1. H2/STO-3G CISD ground energy exactly equals FCI (-1.1372759436 Hartree):
 *      with only 2 electrons, every excitation out of  reference is at most a
 *      double, so CISD and FCI subspaces coincide exactly - a correctness check
 *      for singles/doubles basis-selection logic itself
 *   2. LiH/STO-3G frozen-core CISD (1 frozen + 2 active electrons) exactly
 *      equals frozen-core FCI (-7.882167498160469 Hartree), for same reason: 2
 *      active electrons means CISD=FCI
 *   3. LiH/STO-3G all-electron (4-electron) CISD is a truncated-CI
 *      approximation: matches ci.CISD() reference (-7.882381633123629 Hartree)
 *      and satisfies variational ordering E_RHF >= E_CISD >= E_FCI
 *   4. CISD subspace dimension matches combinatorial formula
 *      1 + n_occ*n_virt + C(n_occ,2)*C(n_virt,2)
 *   5. Invalid-input handling
 */

#include "../physics/cisd.h"
#include "../physics/fci.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "matrix.h"
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

static int binomial(int n, int k) {
  if (k < 0 || k > n) {
    return 0;
  }
  if (k > n - k) {
    k = n - k;
  }

  long long result = 1;
  for (int i = 0; i < k; i++) {
    result = result * (n - i) / (i + 1);
  }

  return (int)result;
}

static void test_h2_cisd_equals_fci(void) {
  printf("Test: H2/STO-3G CISD exactly equals FCI (2 electrons -> every "
         "excitation is at most double)\n");

  double R = 1.4;
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};
  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};
  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);

  molecular_hf_result_t *hf = molecular_rhf(basis, 2, mol, 2, 1e-12, 200);
  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 2, mol);
  double *eri_ao = molecular_eri_tensor(basis, 2);
  double *h_mo = malloc(4 * sizeof(double));
  double *eri_mo = malloc(16 * sizeof(double));

  molecular_ao_to_mo(h_ao, eri_ao, hf->C, 2, h_mo, eri_mo);
  double Enuc = molecule_nuclear_repulsion(mol);

  cisd_result_t *cisd = cisd_solve(2, h_mo, eri_mo, Enuc, 2);
  fci_result_t *fci = fci_solve(2, h_mo, eri_mo, Enuc, 2);

  check(cisd != NULL, "cisd_solve should succeed");
  check(fci != NULL, "fci_solve should succeed");
  if (cisd && fci) {
    check(cisd->dim == fci->dim,
          "CISD and FCI subspace dimensions match exactly for a "
          "2-electron reference (dim = 1 + n_occ * n_virt + C(n_occ,2) * "
          "C(n_virt,2) = C(n_modes,n_electrons))");
    check_close(cisd->ground_energy, -1.1372759436, 1e-6,
                "H2/STO-3G CISD ground energy matches FCI reference");
    check_close(cisd->ground_energy, fci->ground_energy, 1e-9,
                "CISD ground energy exactly equals FCI ground energy");
    check_close(cisd->correlation_energy,
                cisd->ground_energy - cisd->reference_energy, 1e-12,
                "correlation_energy = ground_energy - reference_energy");
    check(cisd->reference_energy > cisd->ground_energy - 1e-9,
          "reference (HF) energy is a worse (higher) variational bound than "
          "CISD ground energy");
  }

  if (cisd) {
    cisd_result_free(cisd);
  }
  if (fci) {
    fci_result_free(fci);
  }
  molecular_hf_result_free(hf);
  cmatrix_free(h_ao);
  free(eri_ao);
  free(h_mo);
  free(eri_mo);
  basis_function_free(h0);
  basis_function_free(h1);
  molecule_free(mol);
}

static void test_lih_setup(basis_function_t **basis_out, molecule_t **mol_out,
                           molecular_hf_result_t **hf_out, double **h_mo_out,
                           double **eri_mo_out, double *enuc_out) {
  double R = 3.015;
  static basis_function_t *li_orbs[5];
  double c_li[3] = {0, 0, 0};
  double c_h[3] = {0, 0, R};

  molint_basis_sto3g_li(c_li, li_orbs);
  basis_function_t *h_orb = molint_basis_sto3g_h(c_h);

  basis_out[0] = li_orbs[0];
  basis_out[1] = li_orbs[1];
  basis_out[2] = li_orbs[2];
  basis_out[3] = li_orbs[3];
  basis_out[4] = li_orbs[4];
  basis_out[5] = h_orb;

  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  *mol_out = molecule_alloc(2, charges, centers);

  *hf_out = molecular_rhf(basis_out, 6, *mol_out, 4, 1e-12, 200);
  cmatrix_t *h_ao = molecular_core_hamiltonian(basis_out, 6, *mol_out);
  double *eri_ao = molecular_eri_tensor(basis_out, 6);
  *h_mo_out = malloc(36 * sizeof(double));
  *eri_mo_out = malloc(6 * 6 * 6 * 6 * sizeof(double));

  molecular_ao_to_mo(h_ao, eri_ao, (*hf_out)->C, 6, *h_mo_out, *eri_mo_out);
  *enuc_out = molecule_nuclear_repulsion(*mol_out);

  cmatrix_free(h_ao);
  free(eri_ao);
}

static void test_lih_frozen_core_cisd_equals_frozen_core_fci(void) {
  printf("Test: LiH/STO-3G frozen-core CISD (2 active electrons) exactly "
         "equals frozen-core FCI\n");

  basis_function_t *basis[6];
  molecule_t *mol;
  molecular_hf_result_t *hf;
  double *h_mo, *eri_mo, enuc;

  test_lih_setup(basis, &mol, &hf, &h_mo, &eri_mo, &enuc);

  cisd_result_t *cisd = cisd_solve_frozen_core(6, 1, h_mo, eri_mo, enuc, 2);
  fci_result_t *fci = fci_solve_frozen_core(6, 1, h_mo, eri_mo, enuc, 2);

  check(cisd != NULL, "cisd_solve_frozen_core should succeed");
  check(fci != NULL, "fci_solve_frozen_core should succeed");
  if (cisd && fci) {
    check(cisd->dim == fci->dim,
          "CISD and FCI subspace dimensions match exactly for a "
          "2-active-electron reference");
    check_close(cisd->ground_energy, -7.882167498160469, 1e-6,
                "frozen-core CISD matches frozen-core FCI reference");
    check_close(cisd->ground_energy, fci->ground_energy, 1e-9,
                "frozen-core CISD ground energy exactly equals frozen-core "
                "FCI ground energy");
  }

  if (cisd) {
    cisd_result_free(cisd);
  }
  if (fci) {
    fci_result_free(fci);
  }
  molecular_hf_result_free(hf);
  free(h_mo);
  free(eri_mo);
  for (int i = 0; i < 5; i++) {
    basis_function_free(basis[i]);
  }
  basis_function_free(basis[5]);
  molecule_free(mol);
}

static void test_lih_all_electron_cisd_variational_ordering(void) {
  printf("Test: LiH/STO-3G all-electron (4-electron) CISD matches reference "
         "and lies strictly between RHF and FCI\n");

  basis_function_t *basis[6];
  molecule_t *mol;
  molecular_hf_result_t *hf;
  double *h_mo, *eri_mo, enuc;

  test_lih_setup(basis, &mol, &hf, &h_mo, &eri_mo, &enuc);

  cisd_result_t *cisd = cisd_solve(6, h_mo, eri_mo, enuc, 4);

  check(cisd != NULL, "cisd_solve should succeed");
  if (cisd) {
    int n_occ = 4, n_virt = 2 * 6 - 4;
    int expected_dim =
        1 + n_occ * n_virt + binomial(n_occ, 2) * binomial(n_virt, 2);

    printf("  CISD ground energy: %.10f Hartree (dim=%d)\n",
           cisd->ground_energy, cisd->dim);

    check(cisd->dim == expected_dim,
          "CISD subspace dimension matches combinatorial formula");
    check_close(cisd->ground_energy, -7.882381633123629, 1e-6,
                "LiH/STO-3G all-electron CISD matches CISD reference");
    check_close(cisd->reference_energy, hf->total_energy, 1e-8,
                "CISD reference-determinant energy matches converged RHF "
                "total energy");

    check(cisd->reference_energy > cisd->ground_energy + 1e-6,
          "CISD ground energy is strictly below (better than) RHF "
          "reference energy");
    check(cisd->ground_energy > -7.8823949575 + 1e-6,
          "CISD ground energy is strictly above (a worse variational bound "
          "than) full-active-space FCI ground energy, since CISD's subspace is "
          "a strict subset of FCI's");

    cisd_result_free(cisd);
  }

  molecular_hf_result_free(hf);
  free(h_mo);
  free(eri_mo);
  for (int i = 0; i < 5; i++) {
    basis_function_free(basis[i]);
  }
  basis_function_free(basis[5]);
  molecule_free(mol);
}

static void test_invalid_inputs_rejected(void) {
  printf("Test: invalid inputs are rejected cleanly\n");

  const double h_mo[4] = {0};
  const double eri_mo[16] = {0};

  check(cisd_solve(0, h_mo, eri_mo, 0.0, 2) == NULL,
        "n_spatial<1 should be rejected");
  check(cisd_solve(2, NULL, eri_mo, 0.0, 2) == NULL,
        "NULL h_mo should be rejected");
  check(cisd_solve(2, h_mo, NULL, 0.0, 2) == NULL,
        "NULL eri_mo should be rejected");
  check(cisd_solve(2, h_mo, eri_mo, 0.0, -1) == NULL,
        "negative n_electrons should be rejected");
  check(cisd_solve(2, h_mo, eri_mo, 0.0, 5) == NULL,
        "n_electrons > 2*n_spatial should be rejected");

  check(cisd_solve_frozen_core(2, 2, h_mo, eri_mo, 0.0, 0) == NULL,
        "n_frozen >= n_spatial should be rejected");
  check(cisd_solve_frozen_core(2, -1, h_mo, eri_mo, 0.0, 0) == NULL,
        "negative n_frozen should be rejected");
}

int main(void) {
  test_h2_cisd_equals_fci();
  test_lih_frozen_core_cisd_equals_frozen_core_fci();
  test_lih_all_electron_cisd_variational_ordering();
  test_invalid_inputs_rejected();

  if (failures == 0) {
    printf("\nAll test_cisd checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
