/*
 * Test: Variational Quantum Eigensolver (hardware-efficient ansatz +
 * coordinate-descent optimizer).
 *
 * 1. Single-qubit H = a*X + b*Z: exact ground energy is -\sqrt{a^2 + b^2}
 *    (2x2 Hermitian, closed form).
 * 2. Transverse-field Ising model (n_qubits=3): cross-checked against exact
 *    diagonalization via cmatrix_eigh_complex, since TFIM's ground energy has
 *    no simple closed form for a small open chain.
 * 3. vqe_expectation / vqe_energy fixture checks at fixed (untrained)
 *    parameters.
 * 4. Invalid-input handling.
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../core/sparse.h"
#include "../physics/spin_chain.h"
#include "../physics/vqe.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static int hermitian_ok(const cmatrix_t *H) {
  for (int i = 0; i < H->nrows; i++) {
    for (int j = 0; j < H->ncols; j++) {
      complex_t diff = c_sub(CMAT(H, i, j), c_conj(CMAT(H, j, i)));
      if (sqrt(c_abs2(diff)) > 1e-10) {
        return 0;
      }
    }
  }

  return 1;
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

  // \theta = \pi flips |0> -> |1> (up to phase) via RY(pi), so <1|H|1>
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

// Independent (test-only) reference: builds a Pauli-string Hamiltonian via
// explicit Kronecker products
static void kron_into(const complex_t *a, int da, const complex_t *b, int db,
                      complex_t *out) {
  int dout = da * db;
  for (int i = 0; i < da; i++) {
    for (int j = 0; j < da; j++) {
      complex_t aij = a[i * da + j];
      for (int p = 0; p < db; p++) {
        for (int q = 0; q < db; q++) {
          out[(i * db + p) * dout + (j * db + q)] = c_mul(aij, b[p * db + q]);
        }
      }
    }
  }
}

static cmatrix_t *reference_pauli_hamiltonian(int n_qubits,
                                              const char *const *strings,
                                              const double *coeffs,
                                              int n_terms) {
  complex_t I2[4] = {c_new(1, 0), c_zero(), c_zero(), c_new(1, 0)};
  complex_t X2[4] = {c_zero(), c_new(1, 0), c_new(1, 0), c_zero()};
  complex_t Y2[4] = {c_zero(), c_new(0, -1), c_new(0, 1), c_zero()};
  complex_t Z2[4] = {c_new(1, 0), c_zero(), c_zero(), c_new(-1, 0)};

  int dim = 1 << n_qubits;
  cmatrix_t *H = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    H->data[i] = c_zero();
  }

  complex_t *acc = malloc((size_t)dim * dim * sizeof *acc);
  complex_t *tmp = malloc((size_t)dim * dim * sizeof *tmp);

  for (int t = 0; t < n_terms; t++) {
    int cur_dim = 1;
    acc[0] = c_new(1.0, 0.0);

    for (int k = 0; k < n_qubits; k++) {
      const complex_t *op;
      switch (strings[t][k]) {
      case 'I':
        op = I2;
        break;
      case 'X':
        op = X2;
        break;
      case 'Y':
        op = Y2;
        break;
      default:
        op = Z2;
        break;
      }

      kron_into(acc, cur_dim, op, 2, tmp);
      cur_dim *= 2;
      memcpy(acc, tmp, (size_t)cur_dim * cur_dim * sizeof *acc);
    }

    for (int idx = 0; idx < dim * dim; idx++) {
      H->data[idx] = c_add(H->data[idx], c_scale(acc[idx], coeffs[t]));
    }
  }

  free(acc);
  free(tmp);

  return H;
}

static void test_pauli_hamiltonian(void) {
  printf("test_pauli_hamiltonian:\n");

  // Single-term sanity checks against hand-known 2x2/4x4 matrices.
  {
    const char *s = "X";
    double c = 1.0;
    cmatrix_t *H = vqe_build_pauli_hamiltonian(1, &s, &c, 1);
    check_true(H != NULL, "single-qubit X builds");
    if (H) {
      check_close(CMAT(H, 0, 1).re, 1.0, 1e-12, "X: <0|X|1>=1");
      check_close(CMAT(H, 1, 0).re, 1.0, 1e-12, "X: <1|X|0>=1");
      check_close(CMAT(H, 0, 0).re, 0.0, 1e-12, "X: <0|X|0>=0");

      cmatrix_free(H);
    }
  }
  {
    const char *s = "Y";
    double c = 1.0;
    cmatrix_t *H = vqe_build_pauli_hamiltonian(1, &s, &c, 1);
    check_true(H != NULL, "single-qubit Y builds");
    if (H) {
      check_close(CMAT(H, 0, 1).im, -1.0, 1e-12, "Y: <0|Y|1>=-i");
      check_close(CMAT(H, 1, 0).im, 1.0, 1e-12, "Y: <1|Y|0>=+i");

      cmatrix_free(H);
    }
  }
  {
    const char *s = "Z";
    double c = 1.0;
    cmatrix_t *H = vqe_build_pauli_hamiltonian(1, &s, &c, 1);
    check_true(H != NULL, "single-qubit Z builds");
    if (H) {
      check_close(CMAT(H, 0, 0).re, 1.0, 1e-12, "Z: <0|Z|0>=+1");
      check_close(CMAT(H, 1, 1).re, -1.0, 1e-12, "Z: <1|Z|1>=-1");

      cmatrix_free(H);
    }
  }

  // Multi-term, multi-qubit: cross-check against the independent
  // Kronecker-product reference for several random-ish Pauli sums.
  {
    const char *strings[] = {"XI", "IY", "ZZ", "XY"};
    const double coeffs[] = {0.7, -0.3, 1.1, 0.4};
    cmatrix_t *H = vqe_build_pauli_hamiltonian(2, strings, coeffs, 4);
    cmatrix_t *Href = reference_pauli_hamiltonian(2, strings, coeffs, 4);

    double max_err = 0.0;
    for (int i = 0; i < 4 * 4; i++) {
      double e = sqrt(c_abs2(c_sub(H->data[i], Href->data[i])));
      if (e > max_err) {
        max_err = e;
      }
    }
    check_close(max_err, 0.0, 1e-12,
                "2-qubit multi-term sum matches independent Kronecker "
                "reference");
    check_true(hermitian_ok(H), "2-qubit multi-term sum is Hermitian");

    cmatrix_free(H);
    cmatrix_free(Href);
  }

  // TFIM cross-check: build the same model both via vqe_build_tfim
  // (hand-rolled bit manipulation) and via this general Pauli-string
  // layer, and confirm they agree element-by-element.
  {
    int n = 4;
    double J = 1.0, h = 0.5;
    cmatrix_t *H_bits = vqe_build_tfim(n, J, h);

    const char *strings[7];
    double coeffs[7];
    char buf[7][5];
    int nt = 0;
    for (int k = 0; k < n - 1; k++) {
      for (int q = 0; q < n; q++) {
        buf[nt][q] = (q == k || q == k + 1) ? 'Z' : 'I';
      }

      buf[nt][n] = '\0';
      strings[nt] = buf[nt];
      coeffs[nt] = -J;
      nt++;
    }
    for (int k = 0; k < n; k++) {
      for (int q = 0; q < n; q++) {
        buf[nt][q] = (q == k) ? 'X' : 'I';
      }

      buf[nt][n] = '\0';
      strings[nt] = buf[nt];
      coeffs[nt] = -h;
      nt++;
    }

    cmatrix_t *H_pauli = vqe_build_pauli_hamiltonian(n, strings, coeffs, nt);

    double max_err = 0.0;
    for (int i = 0; i < (1 << n) * (1 << n); i++) {
      double e = sqrt(c_abs2(c_sub(H_bits->data[i], H_pauli->data[i])));
      if (e > max_err) {
        max_err = e;
      }
    }
    check_close(max_err, 0.0, 1e-12,
                "TFIM (n=4) matches vqe_build_tfim's own bit-manipulation "
                "construction exactly, term-by-term via Pauli strings");

    cmatrix_free(H_bits);
    cmatrix_free(H_pauli);
  }

  // NOTE: XXZ Heisenberg chain:
  //  H = \sum_j [(Jz/4) Z_jZ_{j+1} + (Jxy/4)(X_jX_{j+1} + Y_jY_{j+1})]
  // is spin-1/2 <-> qubit mapping (S^z=Z/2, S+S-+S-S+ = (XX+YY)/2) of the exact
  // same XXZ Hamiltonian spin_sector_hamiltonian() diagonalizes via completely
  // independent machinery (bit-manipulation + sparse matrices + symmetry
  // sectors, not dense tensor products). Comparing their ground energies is a
  // strongc cross-check between two unrelated code paths in this codebase.
  {
    int N = 4;
    double Jz = 1.3, Jxy = 0.7;

    const char *strings[12];
    double coeffs[12];
    char buf[12][5];
    int nt = 0;
    for (int j = 0; j < N; j++) {
      int jn = (j + 1) % N;
      const char ops[3] = {'Z', 'X', 'Y'};
      const double base[3] = {Jz / 4.0, Jxy / 4.0, Jxy / 4.0};
      for (int o = 0; o < 3; o++) {
        for (int q = 0; q < N; q++) {
          buf[nt][q] = (q == j || q == jn) ? ops[o] : 'I';
        }

        buf[nt][N] = '\0';
        strings[nt] = buf[nt];
        coeffs[nt] = base[o];
        nt++;
      }
    }

    cmatrix_t *H = vqe_build_pauli_hamiltonian(N, strings, coeffs, nt);
    check_true(hermitian_ok(H), "XXZ (N=4) Pauli-sum Hamiltonian is Hermitian");

    cmatrix_t *Hc = cmatrix_copy(H);
    eigen_t *eig = cmatrix_eigh_complex(Hc);
    cmatrix_free(Hc);

    double e0_pauli = eig ? eig->eigenvalues[0] : 1e300;
    if (eig) {
      eigen_free(eig);
    }

    // Best (most negative) ground energy across every (nup, k) sector, matching
    // spin_sector_hamiltonian's own symmetry decomposition - full spectrum is
    // the union of every sector's spectrum.
    double e0_sector = 1e300;
    for (int nup = 0; nup <= N; nup++) {
      for (int k = 0; k < N; k++) {
        spin_sector_t *sec = spin_sector_build(N, nup, k);
        if (!sec) {
          continue;
        }

        sparse_matrix_t *Hs = spin_sector_hamiltonian(sec, Jxy, Jz, /*pbc=*/1);

        // Dense diagonalization, not Lanczos: sectors are tiny for N=4 (up to a
        // handful of states), and Lanczos with max_iter == dim can hit early
        // breakdown on such small matrices - dense is exact and cheap here, and
        // this cross-check needs to be authoritative.
        int sdim = sec->dim;
        if (sdim <= 0) {
          sparse_free(Hs);
          spin_sector_free(sec);

          continue;
        }

        cmatrix_t *Hd = cmatrix_alloc(sdim, sdim);
        for (int ii = 0; ii < sdim * sdim; ii++) {
          Hd->data[ii] = c_zero();
        }
        for (int ii = 0; ii < sdim; ii++) {
          for (int idx = Hs->row_ptr[ii]; idx < Hs->row_ptr[ii + 1]; idx++) {
            int jj = Hs->col_ind[idx];

            CMAT(Hd, ii, jj) = c_add(CMAT(Hd, ii, jj), Hs->values[idx]);
          }
        }

        eigen_t *sector_eig = cmatrix_eigh_complex(Hd);
        if (sector_eig && sector_eig->eigenvalues[0] < e0_sector) {
          e0_sector = sector_eig->eigenvalues[0];
        }

        if (sector_eig) {
          eigen_free(sector_eig);
        }
        cmatrix_free(Hd);

        sparse_free(Hs);
        spin_sector_free(sec);
      }
    }

    check_close(
        e0_pauli, e0_sector, 1e-8,
        "XXZ (N=4) ground energy: Pauli-string dense diagonalization matches "
        "spin_chain.c's independent sector-by-sector exact diagonalization");

    cmatrix_free(H);
  }

  // Invalid input handling.
  {
    const char *s1 = "X";
    double c1 = 1.0;
    check_true(vqe_build_pauli_hamiltonian(0, &s1, &c1, 1) == NULL,
               "n_qubits<1 rejected");
    check_true(vqe_build_pauli_hamiltonian(1, &s1, &c1, 0) == NULL,
               "n_terms<1 rejected");
    check_true(vqe_build_pauli_hamiltonian(1, NULL, &c1, 1) == NULL,
               "NULL pauli_strings rejected");
    check_true(vqe_build_pauli_hamiltonian(1, &s1, NULL, 1) == NULL,
               "NULL coefficients rejected");

    const char *bad_char = "Q";
    check_true(vqe_build_pauli_hamiltonian(1, &bad_char, &c1, 1) == NULL,
               "invalid Pauli character rejected");

    const char *too_short = "X";
    check_true(vqe_build_pauli_hamiltonian(2, &too_short, &c1, 1) == NULL,
               "string shorter than n_qubits rejected");

    const char *too_long = "XXX";
    check_true(vqe_build_pauli_hamiltonian(2, &too_long, &c1, 1) == NULL,
               "string longer than n_qubits rejected");
  }
}

int main(void) {
  test_expectation_fixtures();
  test_invalid_input();
  test_vqe_single_qubit();
  test_vqe_tfim();
  test_pauli_hamiltonian();

  if (failures == 0) {
    printf("\nAll test_vqe checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
