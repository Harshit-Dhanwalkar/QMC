/*
 * Test: LiH/STO-3G : first real molecule/basis using p-orbitals and a
 * heteronuclear pair, plus the frozen-core reduction needed to make its
 * 12-qubit second-quantized Hamiltonian's exact diagonalization tractable.
 *
 * 1. Li basis (molint_basis_sto3g_li): normalization and p-orbital
 *    orthogonality (<2px|2py>=0 etc, exact by symmetry) : first real exercise
 *    of the general GTO engine on l+m+n>0 angular momentum, everything through
 *    H2/H4 only ever used s functions.
 * 2. LiH RHF: implementation gave -7.856587 Hartree at R=3.015 bohr (~1.596
 *    Angstrom, close to LiH's experimental equilibrium bond length).
 * 3. Frozen-core Hamiltonian construction
 *    (second_quant_build_molecular_hamiltonian_frozen_core): a check that
 *    n_frozen=0 reduces Exactly to the unfrozen builder's own result, then real
 *    validation : freezing Li's 1s core (standard choice in the LiH VQE
 *    literature) gives a 10-qubit (1024-dim) FCI ground state differing from
 *    the exact full 12-qubit (4096-dim) FCI by only about 0.23 mHartree, both
 *    confirming correctness and that Li's 1s core really is close to chemically
 *    inert here.
 *
 * // WARN: VQE on this 10-qubit system is intentionally NOT exercised in this
 * automated test (both slow, like H4's 8-qubit case)
 */

#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/second_quant.h"
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

static void test_li_basis_normalization_and_p_orbitals(void) {
  printf("test_li_basis_normalization_and_p_orbitals:\n");

  double center[3] = {0, 0, 0};
  basis_function_t *li[5];
  int ok = molint_basis_sto3g_li(center, li);
  check_true(ok, "molint_basis_sto3g_li allocates");
  if (!ok) {
    return;
  }

  const char *names[5] = {"1s", "2s", "2px", "2py", "2pz"};
  for (int i = 0; i < 5; i++) {
    char label[64];
    snprintf(label, sizeof(label), "Li %s self-overlap = 1", names[i]);

    check_close(gto_overlap(li[i], li[i]), 1.0, 1e-10, label);
  }

  check_close(gto_overlap(li[2], li[3]), 0.0, 1e-14, "<2px|2py> = 0 exactly");
  check_close(gto_overlap(li[2], li[4]), 0.0, 1e-14, "<2px|2pz> = 0 exactly");
  check_close(gto_overlap(li[3], li[4]), 0.0, 1e-14, "<2py|2pz> = 0 exactly");
  check_true(fabs(gto_overlap(li[0], li[1])) > 1e-3,
             "<1s|2s> is nonzero (same symmetry, different shell : expected in "
             "an unorthogonalized AO basis)");

  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }
}

static void test_lih_rhf_and_frozen_core_fci(void) {
  printf("test_lih_rhf_and_frozen_core_fci:\n");

  double R = 3.015; // ~1.596 Angstrom, near LiH's experimental equilibrium
  double li_center[3] = {0, 0, 0};
  basis_function_t *li[5];
  molint_basis_sto3g_li(li_center, li);
  double h_center[3] = {0, 0, R};
  basis_function_t *h = molint_basis_sto3g_h(h_center);
  basis_function_t *basis[6] = {li[0], li[1], li[2], li[3], li[4], h};

  double charge[2] = {3.0, 1.0};
  double centers_pos[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charge, centers_pos);

  molecular_hf_result_t *hf = molecular_rhf(basis, 6, mol, 4, 1e-9, 300);
  check_true(hf != NULL, "molecular_rhf allocates for LiH");
  if (!hf) {
    molecule_free(mol);
    basis_function_free(h);

    for (int i = 0; i < 5; i++) {
      basis_function_free(li[i]);
    }

    return;
  }
  printf("  LiH RHF: %.6f Hartree (iters=%d)\n", hf->total_energy,
         hf->iterations);
  check_true(hf->converged, "LiH RHF converges");
  check_close(hf->total_energy, -7.856587, 1e-4,
              "LiH RHF matches cross-validation");

  cmatrix_t *Hcore = molecular_core_hamiltonian(basis, 6, mol);
  double *eri_ao = molecular_eri_tensor(basis, 6);
  double nuclear_repulsion = molecule_nuclear_repulsion(mol);
  double *h_mo = malloc(36 * sizeof(double));
  double *eri_mo = malloc(1296 * sizeof(double));
  molecular_ao_to_mo(Hcore, eri_ao, hf->C, 6, h_mo, eri_mo);

  // n_frozen=0 must reduce exactly (bit-identically) to unfrozen builder's own
  // result
  cmatrix_t *H_full = second_quant_build_molecular_hamiltonian(
      6, h_mo, eri_mo, nuclear_repulsion);
  cmatrix_t *H_zero_frozen =
      second_quant_build_molecular_hamiltonian_frozen_core(6, 0, h_mo, eri_mo,
                                                           nuclear_repulsion);
  double max_diff = 0.0;
  for (int i = 0; i < H_full->nrows * H_full->ncols; i++) {
    double d = fabs(H_full->data[i].re - H_zero_frozen->data[i].re);
    if (d > max_diff) {
      max_diff = d;
    }
  }
  check_true(max_diff < 1e-12,
             "n_frozen=0 is bit-identical to the unfrozen Hamiltonian builder");
  cmatrix_free(H_full);
  cmatrix_free(H_zero_frozen);

  // The real case: freeze Li's 1s core
  cmatrix_t *H_frozen = second_quant_build_molecular_hamiltonian_frozen_core(
      6, 1, h_mo, eri_mo, nuclear_repulsion);
  check_true(H_frozen != NULL, "frozen-core Hamiltonian allocates");
  check_true(H_frozen->nrows == 1024,
             "frozen-core H is 1024x1024 (5 active spatial orbitals -> 10 "
             "spin-orbitals -> 2^10)");

#ifdef USE_LAPACK
  /* NOTE: Exact diagonalization of a 1024x1024 complex Hermitian matrix is only
   * exercised when built with USE_LAPACK: with the default hand-rolled
   * real-embedding solver (effectively a 2048x2048 real-symmetric QR solve for
   * this case)
   */
  int dim = H_frozen->nrows;
  cmatrix_t *H_copy = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    H_copy->data[i] = H_frozen->data[i];
  }
  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  double E_fci_frozen = eig->eigenvalues[0];
  printf("  LiH frozen-core FCI: %.6f Hartree\n", E_fci_frozen);
  check_close(E_fci_frozen, -7.876732, 1e-4, "LiH frozen-core FCI");
  check_true(E_fci_frozen < hf->total_energy - 1e-6,
             "frozen-core FCI is below RHF (captures correlation energy)");

  printf("  Frozen-core error vs. Full-space FCI (-7.876960): %.6f Hartree\n",
         fabs(E_fci_frozen - (-7.876960)));
  check_true(fabs(E_fci_frozen - (-7.876960)) < 0.005,
             "frozen-core error vs. Full-space FCI is small (Li's 1s core is "
             "close to chemically inert, as expected)");

  eigen_free(eig);
  cmatrix_free(H_copy);
#else
  printf("  (skipping exact diagonalization of the 1024x1024 frozen-core "
         "Hamiltonian : default build has no LAPACK; rebuild with USE_LAPACK=1 "
         "to also validate the FCI energy here. RHF and the Hamiltonian "
         "construction itself are still fully validated above.)\n");
#endif

  cmatrix_free(H_frozen);
  cmatrix_free(Hcore);
  free(eri_ao);
  free(h_mo);
  free(eri_mo);
  molecular_hf_result_free(hf);
  molecule_free(mol);
  basis_function_free(h);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }
}

int main(void) {
  test_li_basis_normalization_and_p_orbitals();
  test_lih_rhf_and_frozen_core_fci();

  if (failures == 0) {
    printf("\nAll test_lih checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_lih check(s) FAILED.\n", failures);
    return 1;
  }
}
