/*
 * Test: UHF (Unrestricted Hartree-Fock) on the Li atom, STO-3G basis.
 *
 * NOTE: Li (Z=3, doublet ground state: 2-\alpha electrons, 1-\beta) is the
 * smallest open-shell system reachable with the existing STO-3G basis builders
 * (molint_basis_sto3g_li), and has no even-electron-count RHF analogue on the
 * same atom to compare against directly.
 *
 * <S^2> = 0.75 exactly is also an independent physical check: it is Sz(Sz+1)
 * for a pure doublet (Sz=1/2) with zero spin contamination, which is expected
 * here since Li's unpaired 2s electron is well-separated in energy from other
 * configurations at this basis-set/geometry.
 */

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
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

static void test_li_atom_uhf(void) {
  printf("test_li_atom_uhf:\n");

  double center[3] = {0, 0, 0};
  basis_function_t *li[5];
  int ok = molint_basis_sto3g_li(center, li);
  check_true(ok, "molint_basis_sto3g_li allocates");
  if (!ok) {
    return;
  }

  const double charge[1] = {3.0};
  double centers_pos[1][3] = {{0, 0, 0}};
  molecule_t *mol = molecule_alloc(1, charge, centers_pos);

  // Doublet ground state: 2-\alpha, 1-\beta
  molecular_uhf_result_t *uhf = molecular_uhf(li, 5, mol, 2, 1, 1e-10, 300);
  check_true(uhf != NULL, "molecular_uhf allocates for Li atom");
  if (!uhf) {
    molecule_free(mol);

    for (int i = 0; i < 5; i++) {
      basis_function_free(li[i]);
    }

    return;
  }

  printf("  Li UHF: %.10f Hartree (iters=%d)\n", uhf->total_energy,
         uhf->iterations);
  check_true(uhf->converged, "Li UHF converges");
  check_close(uhf->total_energy, -7.3155259813, 1e-6,
              "Li UHF matches independent reference");
  check_close(uhf->spin_squared, 0.75, 1e-6,
              "<S^2> = Sz(Sz+1) exactly: no spin contamination for Li");

  // NOTE: n_alpha=n_beta edge case: closed-shell UHF on Li+ (2 electrons,
  // singlet) should reduce to \alpha and \beta MOs being numerically identical
  // (both spin channels see the same Fock operator when densities are equal),
  // and its energy should match closed-shell RHF run on the same 2-electron
  // system as an internal consistency check.
  molecular_uhf_result_t *uhf_liplus =
      molecular_uhf(li, 5, mol, 1, 1, 1e-10, 300);
  molecular_hf_result_t *rhf_liplus = molecular_rhf(li, 5, mol, 2, 1e-10, 300);
  check_true(uhf_liplus != NULL && rhf_liplus != NULL,
             "UHF/RHF both allocate for Li+ (2-electron closed shell)");
  if (uhf_liplus && rhf_liplus) {
    check_close(uhf_liplus->total_energy, rhf_liplus->total_energy, 1e-8,
                "UHF with n_alpha=n_beta matches RHF on the same 2e- system");
    check_close(uhf_liplus->spin_squared, 0.0, 1e-8,
                "<S^2> = 0 exactly for a closed-shell singlet");
  }

  molecular_uhf_result_free(uhf_liplus);
  molecular_hf_result_free(rhf_liplus);
  molecular_uhf_result_free(uhf);
  molecule_free(mol);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }
}

int main(void) {
  test_li_atom_uhf();

  if (failures == 0) {
    printf("\nAll test_uhf checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_uhf check(s) FAILED.\n", failures);
    return 1;
  }
}
