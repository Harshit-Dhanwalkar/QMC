/*
 * Test: noisy VQE (density-matrix ansatz simulation under T1/T2 Lindblad noise,
 * layered on top of vqe.c's hardware-efficient ansatz) and zero-noise
 * extrapolation (ZNE).
 *
 * 1. Zero-noise limit: \gamma1=\gamma2=0 must reproduce vqe.c's pure-state
 *   vqe_energy exactly, since a decoherence rate of zero degenerates the
 *   density-matrix evolution to ordinary unitary conjugation.
 * 2. Physicality: Tr(\rho)=1 and purity in (0,1] at nonzero noise, and purity
 *   strictly decreases as gamma increases (state becomes more mixed).
 * 3. ZNE sanity: on a case with a controlled amount of noise, the extrapolated
 *   energy is closer to the true (noiseless) ground energy than the raw noisy
 *   energy at the base noise level.
 * 4. Invalid-input handling.
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/vqe.h"
#include "../physics/vqe_noisy.h"
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

static double density_trace(const cmatrix_t *rho) {
  double tr = 0.0;
  for (int i = 0; i < rho->nrows; i++) {
    tr += CMAT(rho, i, i).re;
  }

  return tr;
}

static double density_purity(const cmatrix_t *rho) {
  double p = 0.0;
  for (int i = 0; i < rho->nrows * rho->ncols; i++) {
    p += c_abs2(rho->data[i]);
  }

  return p;
}

static void test_zero_noise_limit(void) {
  printf("test_zero_noise_limit:\n");

  int n_qubits = 3, n_layers = 2;
  int n_params = n_qubits * n_layers;
  const double theta[6] = {0.3, -0.7, 1.1, 0.2, -0.4, 0.9};

  cmatrix_t *H = vqe_build_tfim(n_qubits, 1.0, 0.5);

  double e_pure = vqe_energy(n_qubits, n_layers, theta, H);
  double e_noisy_zero =
      vqe_noisy_energy(n_qubits, n_layers, theta, H, 0.0, 0.0, 0.3);

  printf("  pure=%.10f noisy(\\gamma=0)=%.10f\n", e_pure, e_noisy_zero);
  check_close(
      e_noisy_zero, e_pure, 1e-9,
      "vqe_noisy_energy(\\gamma1=\\gamma2=0) matches noiseless vqe_energy");

  cmatrix_t *rho =
      vqe_noisy_prepare_density(n_qubits, n_layers, theta, 0.0, 0.0, 0.3);
  check_true(rho != NULL, "vqe_noisy_prepare_density allocates");
  if (rho) {
    check_close(density_trace(rho), 1.0, 1e-9,
                "Tr(\\rho)=1 at zero noise (still a valid pure state)");
    check_close(density_purity(rho), 1.0, 1e-9,
                "purity=1 at zero noise (still a pure state)");
    cmatrix_free(rho);
  }

  cmatrix_free(H);
  (void)n_params;
}

static void test_physicality_and_purity_decay(void) {
  printf("test_physicality_and_purity_decay:\n");

  int n_qubits = 2, n_layers = 2;
  const double theta[4] = {0.8, -1.2, 0.5, 1.7};

  const double gammas[3] = {0.01, 0.05, 0.15};
  double prev_purity = 1.0 + 1.0; /* sentinel, larger than any valid purity */

  for (int k = 0; k < 3; k++) {
    cmatrix_t *rho = vqe_noisy_prepare_density(n_qubits, n_layers, theta,
                                               gammas[k], gammas[k] * 0.6, 0.4);
    check_true(rho != NULL, "vqe_noisy_prepare_density allocates (k)");
    if (!rho) {
      continue;
    }

    double tr = density_trace(rho);
    double purity = density_purity(rho);
    printf("  gamma1=%.3f: Tr(\\rho)=%.10f purity=%.6f\n", gammas[k], tr,
           purity);

    check_close(tr, 1.0, 1e-6, "Tr(\\rho)=1 (trace-preserving GKSL evolution)");
    check_true(purity <= 1.0 + 1e-9 && purity > 0.0,
               "purity in (0, 1] (valid mixed-state density matrix)");
    check_true(purity < prev_purity,
               "purity strictly decreases as noise (\\gamma) increases");

    prev_purity = purity;
    cmatrix_free(rho);
  }
}

static void test_zne_reduces_error(void) {
  printf("test_zne_reduces_error:\n");

  // Small 2-qubit TFIM: cheap enough for repeated noisy density-matrix
  // simulation across multiple noise scales in a unit test.
  int n_qubits = 2, n_layers = 2;
  cmatrix_t *H = vqe_build_tfim(n_qubits, 1.0, 0.6);

  cmatrix_t *H_copy = cmatrix_copy(H);
  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  cmatrix_free(H_copy);
  check_true(eig != NULL, "exact diagonalization succeeds");
  double E_exact = eig ? eig->eigenvalues[0] : 0.0;

  vqe_result_t r = vqe_run(n_qubits, n_layers, H, 8, M_PI, 4242ULL);
  check_true(r.theta_opt != NULL, "VQE optimization succeeds");

  double gamma1 = 0.008, gamma2 = 0.004, gate_time = 0.15;
  double raw_c1 = 0.0;
  double zne = vqe_noisy_zne_energy(n_qubits, n_layers, r.theta_opt, H, gamma1,
                                    gamma2, gate_time, 3, &raw_c1);

  double err_raw = fabs(raw_c1 - E_exact);
  double err_zne = fabs(zne - E_exact);
  printf("  exact=%.8f  raw(c=1)=%.8f (err=%.2e)  ZNE=%.8f (err=%.2e)\n",
         E_exact, raw_c1, err_raw, zne, err_zne);
  check_true(err_zne < err_raw,
             "ZNE-extrapolated energy is closer to exact than raw noisy "
             "energy at base noise level");

  free(r.theta_opt);
  if (eig) {
    eigen_free(eig);
  }
  cmatrix_free(H);
}

static void test_invalid_input(void) {
  printf("test_invalid_input:\n");

  const double theta[2] = {0.1, 0.2};
  cmatrix_t *H = cmatrix_alloc(2, 2);
  CMAT(H, 0, 0) = c_real(1.0);
  CMAT(H, 0, 1) = c_zero();
  CMAT(H, 1, 0) = c_zero();
  CMAT(H, 1, 1) = c_real(-1.0);

  check_true(vqe_noisy_prepare_density(0, 1, theta, 0.0, 0.0, 0.1) == NULL,
             "n_qubits=0 rejected");
  check_true(vqe_noisy_prepare_density(1, 1, NULL, 0.0, 0.0, 0.1) == NULL,
             "NULL \\theta rejected");
  check_true(vqe_noisy_prepare_density(1, 1, theta, -0.1, 0.0, 0.1) == NULL,
             "negative \\gamma1 rejected");
  check_true(vqe_noisy_prepare_density(1, 1, theta, 0.0, 0.0, -0.1) == NULL,
             "negative gate_time rejected");

  check_close(vqe_noisy_energy(1, 1, theta, NULL, 0.0, 0.0, 0.1), 0.0, 0.0,
              "NULL Hamiltonian returns 0.0");
  check_close(vqe_noisy_energy(2, 1, theta, H, 0.0, 0.0, 0.1), 0.0, 0.0,
              "n_qubits/H dimension mismatch returns 0.0");

  double raw = 0.0;
  check_close(vqe_noisy_zne_energy(1, 1, theta, H, 0.0, 0.0, 0.1, 1, &raw), 0.0,
              0.0, "n_scales<2 returns 0.0");

  cmatrix_free(H);
}

int main(void) {
  test_zero_noise_limit();
  test_physicality_and_purity_decay();
  test_zne_reduces_error();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_vqe_noisy checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
