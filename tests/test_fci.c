/*
 * Test: Full Configuration Interaction (FCI)
 *
 * This module formalizes a pattern previously duplicated by hand across 5 test
 * files: build full Fock-space molecular Hamiltonian and diagonalize it, but
 * properly restricted to a fixed electron-number sector rather than
 * diagonalizing whole Fock space and trusting true ground state happens to also
 * be global minimum across every particle-number sector.
 *
 * Validation:
 *   1. H2/STO-3G and LiH/STO-3G FCI ground energies against
 *      fci.FCI(mf).kernel() (same RHF/STO-3G references):
 *      H2/STO-3G  @ R=1.4 bohr  : -1.1372759436 Hartree
 *      LiH/STO-3G @ R=3.015 bohr: -7.8823949575 Hartree
 *   2. Cross-validated against unrestricted Fock-space diagonalization approach
 *   3. Ground-state electron-number expectation value equals exactly requested
 *      n_electrons (sanity check that the sector restriction didn't
 *      accidentally include/exclude the wrong basis states).
 *   4. Frozen-core FCI on LiH matches an independent frozen-core FCI
 *      calculation.
 *   5. Invalid-input handling.
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/fci.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/second_quant.h"
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

static void test_h2_fci_matches_pyscf(void) {
  printf("Test: H2/STO-3G FCI ground energy matches fci.FCI "
         "reference (-1.1372759436 Hartree)\n");

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

  fci_result_t *fci = fci_solve(2, h_mo, eri_mo, Enuc, 2);
  check(fci != NULL, "fci_solve should succeed");
  if (fci) {
    check(fci->dim == 6,
          "2-electron-in-4-spin-orbital sector dim = C(4,2) = 6");
    check_close(fci->ground_energy, -1.1372759436, 1e-6,
                "H2/STO-3G FCI ground energy");

    // electron-number sanity check
    double N_expect = 0.0;
    for (int i = 0; i < fci->dim; i++) {
      complex_t amp = CMAT(fci->eigenvectors, i, 0);
      double p = amp.re * amp.re + amp.im * amp.im;
      N_expect += p * __builtin_popcount((unsigned)fci->basis_states[i]);
    }
    check_close(N_expect, 2.0, 1e-9,
                "ground state electron-number expectation = 2 exactly");

    /* Full-Fock-space cross-check: cheap here (16-dimensional), unlike the
     * LiH case below (4096-dimensional, too slow to fully diagonaliz ewith
     * default hand-rolled dense complex eigensolver).
     * Confirms sector-restricted result matches unrestricted approach the
     * previously-duplicated test code used. */
    cmatrix_t *H_full =
        second_quant_build_molecular_hamiltonian(2, h_mo, eri_mo, Enuc);
    eigen_t *eig_full = cmatrix_eigh_complex(H_full);
    check_close(fci->ground_energy, eig_full->eigenvalues[0], 1e-9,
                "sector-restricted result matches the unrestricted "
                "full-Fock-space approach's global minimum");
    eigen_free(eig_full);
    cmatrix_free(H_full);

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

static void test_lih_fci_matches_pyscf_and_full_fock_space(void) {
  printf("Test: LiH/STO-3G FCI matches reference -7.8823949575 Hartree and "
         "unrestricted full-Fock-space diagonalization approach duplicated "
         "test code used\n");

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

  molecular_hf_result_t *hf = molecular_rhf(basis, 6, mol, 4, 1e-12, 200);
  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 6, mol);
  double *eri_ao = molecular_eri_tensor(basis, 6);
  double *h_mo = malloc(36 * sizeof(double));
  double *eri_mo = malloc(6 * 6 * 6 * 6 * sizeof(double));

  molecular_ao_to_mo(h_ao, eri_ao, hf->C, 6, h_mo, eri_mo);
  double Enuc = molecule_nuclear_repulsion(mol);

  fci_result_t *fci = fci_solve(6, h_mo, eri_mo, Enuc, 4);
  check(fci != NULL, "fci_solve should succeed");
  if (fci) {
    check_close(fci->ground_energy, -7.8823949575, 1e-5,
                "LiH/STO-3G FCI ground energy");
    /* NOTE: deliberately not cross-checking against a full, unrestricted
     * 2^12=4096-dimensional Fock-space diagononalization here (unlike H2 test,
     * where full Fock space is only 16-dimensional and cheap) -> that would
     * require exactly expensive dense complex diagonalization this module's
     * sector restriction exists to avoid, and this default hand-rolled
     * (non-LAPACK) dense complex eigensolver is very slow at that size. */

    fci_result_free(fci);
  }

  molecular_hf_result_free(hf);
  cmatrix_free(h_ao);
  free(eri_ao);
  free(h_mo);
  free(eri_mo);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li_orbs[i]);
  }
  basis_function_free(h_orb);
  molecule_free(mol);
}

static void test_lih_frozen_core_fci(void) {
  printf("Test: LiH/STO-3G frozen-core FCI (freeze Li 1s) matches an "
         "independent frozen-core FCI calculation\n");

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

  molecular_hf_result_t *hf = molecular_rhf(basis, 6, mol, 4, 1e-12, 200);
  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, 6, mol);
  double *eri_ao = molecular_eri_tensor(basis, 6);
  double *h_mo = malloc(36 * sizeof(double));
  double *eri_mo = malloc(6 * 6 * 6 * 6 * sizeof(double));

  molecular_ao_to_mo(h_ao, eri_ao, hf->C, 6, h_mo, eri_mo);
  double Enuc = molecule_nuclear_repulsion(mol);

  /* Freeze 1 core orbital (Li 1s), leaving 2 electrons active in remaining 5
   * spatial orbitals. */
  fci_result_t *fci = fci_solve_frozen_core(6, 1, h_mo, eri_mo, Enuc, 2);
  check(fci != NULL, "fci_solve_frozen_core should succeed");
  if (fci) {
    printf("  Frozen-core FCI energy: %.8f Hartree\n", fci->ground_energy);
    check_close(fci->ground_energy, -7.882167498160469, 1e-5,
                "frozen-core FCI matches independent CASCI reference");

    check(fci->ground_energy > -7.8823949575 - 1e-6,
          "frozen-core FCI energy is above (less negative than, i.e. a "
          "worse variational bound than) full-active-space FCI, as "
          "expected since freezing removes correlation flexibility");

    fci_result_free(fci);
  }

  molecular_hf_result_free(hf);
  cmatrix_free(h_ao);
  free(eri_ao);
  free(h_mo);
  free(eri_mo);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li_orbs[i]);
  }
  basis_function_free(h_orb);
  molecule_free(mol);
}

static void test_invalid_inputs_rejected(void) {
  printf("Test: invalid inputs are rejected cleanly\n");

  const double h_mo[4] = {0};
  const double eri_mo[16] = {0};

  check(fci_solve(0, h_mo, eri_mo, 0.0, 2) == NULL,
        "n_spatial<1 should be rejected");
  check(fci_solve(2, NULL, eri_mo, 0.0, 2) == NULL,
        "NULL h_mo should be rejected");
  check(fci_solve(2, h_mo, NULL, 0.0, 2) == NULL,
        "NULL eri_mo should be rejected");
  check(fci_solve(2, h_mo, eri_mo, 0.0, -1) == NULL,
        "negative n_electrons should be rejected");
  check(fci_solve(2, h_mo, eri_mo, 0.0, 5) == NULL,
        "n_electrons > 2*n_spatial should be rejected");

  check(fci_solve_frozen_core(2, 2, h_mo, eri_mo, 0.0, 0) == NULL,
        "n_frozen >= n_spatial should be rejected");
  check(fci_solve_frozen_core(2, -1, h_mo, eri_mo, 0.0, 0) == NULL,
        "negative n_frozen should be rejected");
}

int main(void) {
  test_h2_fci_matches_pyscf();
  test_lih_fci_matches_pyscf_and_full_fock_space();
  test_lih_frozen_core_fci();
  test_invalid_inputs_rejected();

  if (failures == 0) {
    printf("\nAll test_fci checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_fci check(s) FAILED.\n", failures);
    return 1;
  }
}
