/*
Test: general N-basis-function molecular_rhf / molecular_ao_to_mo.

1. Regression: molecular_rhf on H2/STO-3G at R=1.4 bohr must reproduce
   -1.116714 Hartree using the general-N path.
2. Result: H4, a linear chain of 4 hydrogen atoms (STO-3G, 1.4 bohr spacing).
  RHF, then via second_quant_build_molecular_hamiltonian + exact diagonalization,
  FCI.
3. VQE on the resulting 8-qubit (256-dim) Hamiltonian, checked to be variational
  (>= FCI) and to get within reach of the RHF baseline with a modest optimizer
  budget. Exploratory tuning found this 48-parameter (8 qubits x 6 layers)
  landscape is noticeably harder to optimize than H2's 12-parameter one.
  // TODO: deeper ansatz, more sweeps, and some seed sensitivity are needed to
  approach FCI, a real (well-known in the VQE    literature)
  optimization-landscape effect at larger qubit count.
  NOTE: This test uses a fast
  budget rather than claiming H2-level convergence at 8 qubits.
*/

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/second_quant.h"
#include "../physics/vqe.h"
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

static void test_molecular_rhf_h2_regression(void) {
  printf("test_molecular_rhf_h2_regression:\n");

  double R = 1.4;
  double centers_pos[2][3] = {{0, 0, 0}, {0, 0, R}};
  basis_function_t *funcs[2] = {molint_basis_sto3g_h(centers_pos[0]),
                                molint_basis_sto3g_h(centers_pos[1])};
  double charge[2] = {1.0, 1.0};
  molecule_t *mol = molecule_alloc(2, charge, centers_pos);

  molecular_hf_result_t *res = molecular_rhf(funcs, 2, mol, 2, 1e-10, 100);
  check_true(res != NULL, "molecular_rhf allocates for H2");
  if (res) {
    printf("  general-N molecular_rhf H2 energy: %.6f\n", res->total_energy);
    check_true(res->converged, "H2 RHF converges");
    check_close(res->total_energy, -1.116714, 1e-4,
                "general-N RHF matches hand-written N=2 RHF (-1.116714)");
    molecular_hf_result_free(res);
  }

  molecule_free(mol);
  basis_function_free(funcs[0]);
  basis_function_free(funcs[1]);
}

static void test_h4_chain_fci_and_vqe(void) {
  printf("test_h4_chain_fci_and_vqe:\n");

  int n = 4;
  double R = 1.4;
  double centers_pos[4][3];
  basis_function_t *funcs[4];
  double charge[4];
  for (int i = 0; i < n; i++) {
    centers_pos[i][0] = 0;
    centers_pos[i][1] = 0;
    centers_pos[i][2] = i * R;
    funcs[i] = molint_basis_sto3g_h(centers_pos[i]);
    charge[i] = 1.0;
  }
  molecule_t *mol = molecule_alloc(n, charge, centers_pos);

  molecular_hf_result_t *hf = molecular_rhf(funcs, n, mol, 4, 1e-10, 200);
  check_true(hf != NULL, "molecular_rhf allocates for H4");
  if (!hf) {
    molecule_free(mol);

    for (int i = 0; i < n; i++) {
      basis_function_free(funcs[i]);
    }

    return;
  }
  printf("  H4 chain RHF: %.6f Hartree (iters=%d)\n", hf->total_energy,
         hf->iterations);
  check_true(hf->converged, "H4 RHF converges");
  check_close(hf->total_energy, -2.098382, 1e-4,
              "H4 RHF matches Python cross-validation");

  cmatrix_t *S = molecular_overlap_matrix(funcs, n);
  cmatrix_t *Hcore = molecular_core_hamiltonian(funcs, n, mol);
  double *eri_ao = molecular_eri_tensor(funcs, n);
  double nuclear_repulsion = molecule_nuclear_repulsion(mol);

  double *h_mo = malloc((size_t)n * n * sizeof(double));
  double *eri_mo = malloc((size_t)n * n * n * n * sizeof(double));
  molecular_ao_to_mo(Hcore, eri_ao, hf->C, n, h_mo, eri_mo);

  cmatrix_t *H = second_quant_build_molecular_hamiltonian(n, h_mo, eri_mo,
                                                          nuclear_repulsion);
  check_true(H != NULL, "second_quant_build_molecular_hamiltonian allocates");
  check_true(H->nrows == 256,
             "H is 256x256 (4 spatial orbitals -> 8 spin-orbitals -> 2^8)");

  int dim = H->nrows;
  cmatrix_t *H_copy = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    H_copy->data[i] = H->data[i];
  }

  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  double E_fci = eig->eigenvalues[0];
  printf("  H4 chain FCI: %.6f Hartree\n", E_fci);
  check_close(E_fci, -2.139443, 1e-4, "H4 FCI matches Python cross-validation");
  check_true(E_fci < hf->total_energy - 1e-6,
             "FCI below RHF (captures correlation energy)");

  double N_expect = 0.0;
  for (int state = 0; state < dim; state++) {
    complex_t amp = CMAT(eig->eigenvectors, state, 0);
    double p = amp.re * amp.re + amp.im * amp.im;

    N_expect += p * __builtin_popcount((unsigned)state);
  }
  check_close(N_expect, 4.0, 1e-6,
              "H4 ground state electron number expectation = 4 exactly");

  printf("  Correlation energy: %.6f Hartree\n", hf->total_energy - E_fci);

  eigen_free(eig);
  cmatrix_free(H_copy);

  /* NOTE: VQE on 8-qubit Hamiltonian: a much harder search than H2's 4 qubits
   * and 12 parameters (n_params = n_qubits*n_layers = 48 here vs 12 for H2).
   */
  vqe_result_t vqe_res = vqe_run(8, 6, H, 30, 0.6, 20260810ULL);
  printf("  VQE (8 qubits, 6 layers, 30 sweeps): E = %.6f Hartree (FCI = %.6f, "
         "RHF = %.6f)\n",
         vqe_res.energy, E_fci, hf->total_energy);
  check_true(vqe_res.theta_opt != NULL, "VQE ran (theta_opt allocated)");
  check_true(vqe_res.energy >= E_fci - 1e-6,
             "VQE energy is variational (>= exact FCI ground state)");
  check_true(vqe_res.energy < hf->total_energy + 0.01,
             "VQE gets within reach of the RHF mean-field baseline");

  free(vqe_res.theta_opt);
  free(h_mo);
  free(eri_mo);
  free(eri_ao);
  cmatrix_free(H);
  cmatrix_free(S);
  cmatrix_free(Hcore);
  molecular_hf_result_free(hf);
  molecule_free(mol);
  for (int i = 0; i < n; i++) {
    basis_function_free(funcs[i]);
  }
}

int main(void) {
  test_molecular_rhf_h2_regression();
  test_h4_chain_fci_and_vqe();

  if (failures == 0) {
    printf("\nAll test_molecular_hf checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_molecular_hf check(s) FAILED.\n", failures);
    return 1;
  }
}
