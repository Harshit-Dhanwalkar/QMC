/*
Test: Variational Quantum Eigensolver (hardware-efficient ansatz +
coordinate-descent optimizer).

1. Single-qubit H = a*X + b*Z: exact ground energy is -\sqrt(a^2 + b^2)
  (2x2 Hermitian, closed form).
2. Transverse-field Ising model (n_qubits=3): cross-checked against exact
  diagonalization via cmatrix_eigh_complex, since TFIM's ground energy has no
  simple closed form for a small open chain.
3. vqe_expectation / vqe_energy fixture checks at fixed (untrained) parameters.
4. Invalid-input handling.
*/

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/vqe.h"
#include <math.h>
#include <stdint.h>
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

static void test_expectation_fixtures(void) {
  printf("test_expectation_fixtures:\n");

  // \theta=0 for a 1-qubit, 1-layer ansatz: RY(0)|0> = |0>, so <0|H|0> =
  // H[0][0]
  cmatrix_t *H = cmatrix_alloc(2, 2);
  CMAT(H, 0, 0) = c_real(1.5);
  CMAT(H, 0, 1) = c_real(0.3);
  CMAT(H, 1, 0) = c_real(0.3);
  CMAT(H, 1, 1) = c_real(-2.1);

  const double theta[1] = {0.0};
  double e = vqe_energy(1, 1, theta, H);
  check_close(e, 1.5, 1e-9, "vqe_energy at \\theta=0 equals H[0][0]");

  // \theta=\pi flips |0> -> |1> (up to phase) via RY(pi), so <1|H|1>
  const double theta_pi[1] = {M_PI};
  double e_pi = vqe_energy(1, 1, theta_pi, H);
  check_close(e_pi, -2.1, 1e-9, "vqe_energy at \\theta=\\pi equals H[1][1]");

  cmatrix_free(H);
}

static void test_vqe_single_qubit(void) {
  printf("test_vqe_single_qubit:\n");

  double a = 1.3, b = 0.7;
  cmatrix_t *H = cmatrix_alloc(2, 2);
  CMAT(H, 0, 0) = c_real(b);
  CMAT(H, 0, 1) = c_real(a);
  CMAT(H, 1, 0) = c_real(a);
  CMAT(H, 1, 1) = c_real(-b);

  double E_exact = -sqrt(a * a + b * b);

  double best_energy = 1e9;
  for (int trial = 0; trial < 5; trial++) {
    vqe_result_t r = vqe_run(1, 2, H, 6, M_PI, 1000ULL + (uint64_t)trial);
    if (r.energy < best_energy) {
      best_energy = r.energy;
    }

    free(r.theta_opt);
  }

  printf("  best of 5 restarts: E=%.8f  exact=%.8f\n", best_energy, E_exact);
  check_close(best_energy, E_exact, 1e-4,
              "VQE (best of 5 restarts) finds 1-qubit ground state");

  cmatrix_free(H);
}

static void test_vqe_tfim(void) {
  printf("test_vqe_tfim:\n");

  int n = 3;
  cmatrix_t *H = vqe_build_tfim(n, 1.0, 0.5);

  cmatrix_t *H_copy = cmatrix_copy(H);
  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  cmatrix_free(H_copy);
  check_true(eig != NULL, "exact diagonalization succeeds");

  double E_exact = eig ? eig->eigenvalues[0] : 0.0;
  if (eig) {
    eigen_free(eig);
  }

  double best_energy = 1e9;
  for (int trial = 0; trial < 8; trial++) {
    vqe_result_t r = vqe_run(n, 3, H, 8, M_PI, 2000ULL + (uint64_t)trial);
    if (r.energy < best_energy) {
      best_energy = r.energy;
    }

    free(r.theta_opt);
  }

  printf("  best of 8 restarts: E=%.8f  exact=%.8f\n", best_energy, E_exact);
  check_close(best_energy, E_exact, 0.05,
              "VQE (best of 8 restarts) approaches the TFIM ground state");
  check_true(best_energy >= E_exact - 1e-6,
             "VQE energy respects the variational bound (E >= exact)");

  cmatrix_free(H);
}

static void test_invalid_input(void) {
  printf("test_invalid_input:\n");

  cmatrix_t *H = cmatrix_alloc(2, 2);
  CMAT(H, 0, 0) = c_real(1.0);
  CMAT(H, 0, 1) = c_zero();
  CMAT(H, 1, 0) = c_zero();
  CMAT(H, 1, 1) = c_real(-1.0);

  check_true(vqe_prepare_ansatz(0, 1, NULL) == NULL, "n_qubits=0 rejected");
  check_true(vqe_prepare_ansatz(1, 1, NULL) == NULL, "NULL \\theta rejected");

  vqe_result_t r1 = vqe_run(1, 1, NULL, 5, M_PI, 1ULL);
  check_true(r1.theta_opt == NULL, "NULL Hamiltonian rejected");

  vqe_result_t r2 =
      vqe_run(2, 1, H, 5, M_PI, 1ULL); // H is 2x2, needs n_qubits=1
  check_true(r2.theta_opt == NULL,
             "Hamiltonian/n_qubits dimension mismatch rejected");

  vqe_result_t r3 = vqe_run(1, 1, H, 5, -1.0, 1ULL);
  check_true(r3.theta_opt == NULL, "non-positive window rejected");

  double e_bad = vqe_expectation(NULL, H);
  check_close(e_bad, 0.0, 0.0, "vqe_expectation(NULL psi) returns 0.0");

  cmatrix_free(H);
}

int main(void) {
  test_expectation_fixtures();
  test_invalid_input();
  test_vqe_single_qubit();
  test_vqe_tfim();

  if (failures == 0) {
    printf("\nAll test_vqe checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
