/*
 * Test: exact diagonalization of the transverse-field Ising model (TFIM) on a
 * periodic ring, cross-checked against its exact Jordan-Wigner + Bogoliubov
 * ground-state energy.
 *
 * Unlike physics/spin_chain.c's XXZ Heisenberg chain, TFIM has a closed-form
 * solution for any N, making it possible to validate against an independent
 * analytic reference rather than another numerical method.
 * ising_exact_ground_energy_per_site() is a from-scratch implementation of that
 * formula.
 *
 * This test checks it against ising_hamiltonian()'s Lanczos ground energy for a
 * spread of N (even And odd, since the antiperiodic-fermion-BC sector exact
 * formula uses turns out to hold for both parities of N) and h/J (including
 * h=0, the classical Ising limit; J=0, decoupled transverse fields; negative h;
 * and h=J, the critical point).
 *
 * Also checks the Z2 (global spin-flip) parity reduction: lower of the +1/-1
 * sector ground energies must match both unreduced Hamiltonian's ground energy
 * AND the exact formula, and the ground state must always sit in the +1 sector
 * for h,J >= 0.
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/sparse.h"
#include "../physics/ising_chain.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  if (err > tol) {
    printf("  FAIL: %s (got=%.12f expected=%.12f err=%.2e > tol=%.2e)\n", label,
           got, expected, err, tol);
    failures++;
  } else {
    printf("  ok:   %s (err=%.2e)\n", label, err);
  }
}

static void check(int cond, const char *label) {
  if (!cond) {
    printf("  FAIL: %s\n", label);
    failures++;
  } else {
    printf("  ok:   %s\n", label);
  }
}

/* Full (unreduced) Hamiltonian's Lanczos ground energy vs exact
 * Jordan-Wigner/Bogoliubov formula, across N=3..10 (even and odd) and a spread
 * of h/J including the critical point, h=0, J=0, and negative h. */
static void test_exact_vs_full_ed(void) {
  printf("-- Full ED vs exact Jordan-Wigner/Bogoliubov ground energy --\n");

  struct {
    int N;
    double J, h;
  } cases[] = {
      {3, 1.0, 0.5},  {3, 1.0, 1.0}, {4, 1.0, 0.0}, {4, 1.0, 1.0},
      {5, 1.0, 0.3},  {5, 1.0, 1.0}, {6, 1.0, 0.7}, {6, 1.0, 1.0},
      {7, 1.0, 1.3},  {8, 1.0, 1.0}, {8, 1.0, 2.0}, {9, 1.0, 1.0},
      {10, 1.0, 1.0}, {6, 0.0, 1.0}, {6, 1.0, 0.0}, {6, 1.0, -0.7},
  };

  for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
    int N = cases[c].N;
    double J = cases[c].J, h = cases[c].h;

    sparse_matrix_t *H = ising_hamiltonian(N, J, h, /*pbc=*/1);
    check(H != NULL, "ising_hamiltonian built");
    if (!H) {
      continue;
    }

    int dim = 1 << N;
    lanczos_result_t *res = lanczos_eigs(H, 1, dim, 1e-12);
    check(res != NULL, "lanczos_eigs converged");
    if (res) {
      double E0_num = res->values[0] / N;
      double E0_exact = ising_exact_ground_energy_per_site(N, J, h);

      char label[128];
      snprintf(label, sizeof label,
               "N=%d J=%.1f h=%.1f: E0/N matches exact formula", N, J, h);

      check_close(E0_num, E0_exact, 1e-8, label);

      lanczos_free(res);
    }

    sparse_free(H);
  }
}

/* Hermiticity of both the full and Z2-parity-reduced Hamiltonians, checked
 * by converting to dense and comparing to its own conjugate transpose.
 * Independent of any of physics, a pure algebraic sanity check. */
static void test_hermiticity(void) {
  printf("\n-- Hermiticity of full and Z2-reduced Hamiltonians --\n");

  int N = 6;
  double J = 1.0, h = 1.1;

  sparse_matrix_t *H = ising_hamiltonian(N, J, h, 1);
  int dim = 1 << N;
  double max_err = 0.0;

  for (int i = 0; i < dim; i++) {
    for (int idx = H->row_ptr[i]; idx < H->row_ptr[i + 1]; idx++) {
      int j = H->col_ind[idx];
      complex_t Hij = H->values[idx];
      // find H[j][i] by linear scan of row j (dim is small here)
      complex_t Hji = c_zero();
      int found = 0;
      for (int idx2 = H->row_ptr[j]; idx2 < H->row_ptr[j + 1]; idx2++) {
        if (H->col_ind[idx2] == i) {
          Hji = H->values[idx2];
          found = 1;

          break;
        }
      }

      double err =
          found ? sqrt(c_abs2(c_sub(Hij, c_conj(Hji)))) : sqrt(c_abs2(Hij));
      if (err > max_err) {
        max_err = err;
      }
    }
  }

  check_close(max_err, 0.0, 1e-12, "full Hamiltonian is Hermitian");
  sparse_free(H);

  sparse_matrix_t *Hp = ising_z2_hamiltonian(N, J, h, 1, +1);
  int dim_half = 1 << (N - 1);
  double max_err_p = 0.0;

  for (int i = 0; i < dim_half; i++) {
    for (int idx = Hp->row_ptr[i]; idx < Hp->row_ptr[i + 1]; idx++) {
      int j = Hp->col_ind[idx];
      complex_t Hij = Hp->values[idx];
      complex_t Hji = c_zero();
      int found = 0;
      for (int idx2 = Hp->row_ptr[j]; idx2 < Hp->row_ptr[j + 1]; idx2++) {
        if (Hp->col_ind[idx2] == i) {
          Hji = Hp->values[idx2];
          found = 1;

          break;
        }
      }

      double err =
          found ? sqrt(c_abs2(c_sub(Hij, c_conj(Hji)))) : sqrt(c_abs2(Hij));
      if (err > max_err_p) {
        max_err_p = err;
      }
    }
  }

  check_close(max_err_p, 0.0, 1e-12, "Z2 parity=+1 Hamiltonian is Hermitian");
  sparse_free(Hp);
}

/* The Z2-reduced +1/-1 sectors' ground energies, merged by taking lower one,
 * must match both the full (unreduced) Hamiltonian's ground energy and the
 * exact formula and ground state must actually sit in the +1 sector, for h,J >=
 * 0. */
static void test_z2_parity_reduction(void) {
  printf("\n-- Z2 parity reduction matches full ED and the exact formula --\n");

  for (int N = 4; N <= 9; N++) {
    double J = 1.0, h = 1.2;

    sparse_matrix_t *Hp = ising_z2_hamiltonian(N, J, h, 1, +1);
    sparse_matrix_t *Hm = ising_z2_hamiltonian(N, J, h, 1, -1);
    check(Hp != NULL && Hm != NULL, "both parity Hamiltonians built");
    if (!Hp || !Hm) {
      sparse_free(Hp);
      sparse_free(Hm);

      continue;
    }

    int dim_half = 1 << (N - 1);
    lanczos_result_t *rp = lanczos_eigs(Hp, 1, dim_half, 1e-12);
    lanczos_result_t *rm = lanczos_eigs(Hm, 1, dim_half, 1e-12);

    double E0_p = rp->values[0] / N;
    double E0_m = rm->values[0] / N;

    char label[128];
    snprintf(label, sizeof label,
             "N=%d: ground state is in the +1 parity sector (E0(+)=%.6f < "
             "E0(-)=%.6f)",
             N, E0_p, E0_m);
    check(E0_p < E0_m, label);

    double E0_exact = ising_exact_ground_energy_per_site(N, J, h);
    snprintf(label, sizeof label,
             "N=%d: parity=+1 sector's ground energy matches exact formula", N);
    check_close(E0_p, E0_exact, 1e-8, label);

    sparse_matrix_t *Hfull = ising_hamiltonian(N, J, h, 1);
    lanczos_result_t *rfull = lanczos_eigs(Hfull, 1, 1 << N, 1e-12);
    double E0_full = rfull->values[0] / N;
    snprintf(label, sizeof label,
             "N=%d: parity=+1 sector's ground energy matches full ED", N);
    check_close(E0_p, E0_full, 1e-8, label);

    lanczos_free(rp);
    lanczos_free(rm);
    lanczos_free(rfull);
    sparse_free(Hp);
    sparse_free(Hm);
    sparse_free(Hfull);
  }
}

/* Basic input-validation guards. */
static void test_invalid_input(void) {
  printf("\n-- Invalid input handling --\n");

  check(ising_hamiltonian(0, 1.0, 1.0, 1) == NULL,
        "ising_hamiltonian rejects N=0");
  check(ising_hamiltonian(31, 1.0, 1.0, 1) == NULL,
        "ising_hamiltonian rejects N=31 (over the uint32_t bit-packing "
        "limit)");
  check(ising_z2_hamiltonian(6, 1.0, 1.0, 1, 0) == NULL,
        "ising_z2_hamiltonian rejects parity=0");
  check(ising_z2_hamiltonian(6, 1.0, 1.0, 1, 2) == NULL,
        "ising_z2_hamiltonian rejects parity=2");
}

/* Entanglement entropy: exact identities that must hold for ANY valid pure
 * state, independent of the physics (S(L_A)=S(N-L_A), and S=0 for a trivial
 * subsystem), plus physically expected trend of decreasing mid-chain
 * entanglement as h moves deep into the paramagnetic phase away from the
 * critical point h=J. */
static void test_entanglement_entropy(void) {
  printf("\n-- Entanglement entropy: exact identities + physical trend --\n");

  int N = 8;
  double J = 1.0;
  double prev_mid = -1.0;
  double hs[] = {0.5, 1.0, 1.5, 3.0, 10.0};

  for (size_t hi = 0; hi < sizeof(hs) / sizeof(hs[0]); hi++) {
    double h = hs[hi];
    sparse_matrix_t *H = ising_hamiltonian(N, J, h, 1);
    int dim = 1 << N;
    lanczos_result_t *res = lanczos_eigs(H, 1, dim, 1e-12);
    check(res != NULL, "ground state found for entanglement entropy test");
    if (!res) {
      sparse_free(H);

      continue;
    }
    const cmatrix_t *psi = res->vectors;

    for (int L_A = 0; L_A <= N; L_A++) {
      double S = ising_entanglement_entropy(psi, N, L_A);
      double S_comp = ising_entanglement_entropy(psi, N, N - L_A);
      char label[128];
      snprintf(
          label, sizeof label,
          "h=%.1f L_A=%d: S(L_A) == S(N-L_A) (subsystem/complement symmetry)",
          h, L_A);
      check_close(S, S_comp, 1e-9, label);

      if (L_A == 0 || L_A == N) {
        snprintf(label, sizeof label,
                 "h=%.1f L_A=%d: trivial subsystem has zero entropy", h, L_A);
        check_close(S, 0.0, 1e-9, label);
      }

      if (L_A == N / 2) {
        if (hi > 0) {
          char label2[192];
          snprintf(
              label2, sizeof label2,
              "h=%.1f: mid-chain entanglement entropy (%.6f) decreased from "
              "previous h (%.6f), moving deeper into the paramagnetic phase",
              h, S, prev_mid);
          check(S < prev_mid, label2);
        }

        prev_mid = S;
      }
    }

    lanczos_free(res);
    sparse_free(H);
  }
}

/* Quench dynamics: exact, parameter-independent identities that must hold for
 * any hermitian H_f and pure state \psi_i, checked here for the TFIM's
 * ising_hamiltonian_dense + ising_loschmidt_echo. L(0)=1 and the no-quench
 * identity are exact to machine precision; short-time Taylor expansion is
 * checked only in the regime where the O(t^4) correction is itself below the
 * comparison tolerance (t=1e-4, 1e-3), and not at larger t where that
 * correction becomes visible */
static void test_quench_dynamics(void) {
  printf("\n-- Quench dynamics: Loschmidt echo exact identities --\n");

  int N = 8;
  double J = 1.0, h_i = 0.5, h_f = 2.0;
  int dim = 1 << N;

  sparse_matrix_t *Hi_sparse = ising_hamiltonian(N, J, h_i, 1);
  lanczos_result_t *gs_i = lanczos_eigs(Hi_sparse, 1, dim, 1e-12);
  check(gs_i != NULL, "pre-quench ground state found");
  if (!gs_i) {
    sparse_free(Hi_sparse);

    return;
  }

  const cmatrix_t *psi_i = gs_i->vectors;

  cmatrix_t *Hf_dense = ising_hamiltonian_dense(N, J, h_f, 1);
  check(Hf_dense != NULL, "post-quench dense Hamiltonian built");

  eigen_t *eig_f = cmatrix_eigh_complex(Hf_dense);
  check(eig_f != NULL, "post-quench Hamiltonian diagonalized");

  if (eig_f) {
    double L0 = ising_loschmidt_echo(eig_f, psi_i, 0.0);
    check_close(L0, 1.0, 1e-9, "L(0) = 1 exactly");

    /* Var(H_f) in the pre-quench state, computed independently of
     * ising_loschmidt_echo via direct matrix-vector products, for short-time
     * Taylor cross-check. */
    cmatrix_t *Hpsi = cmatrix_alloc(dim, 1);
    for (int i = 0; i < dim; i++) {
      complex_t s = c_zero();

      for (int j = 0; j < dim; j++) {
        s = c_add(s, c_mul(CMAT(Hf_dense, i, j), CMAT(psi_i, j, 0)));
      }

      CMAT(Hpsi, i, 0) = s;
    }

    double E_mean = 0.0, E2_mean = 0.0;
    for (int i = 0; i < dim; i++) {
      complex_t e_i = c_mul(c_conj(CMAT(psi_i, i, 0)), CMAT(Hpsi, i, 0));

      E_mean += e_i.re;
      E2_mean += c_abs2(CMAT(Hpsi, i, 0));
    }

    double var_Hf = E2_mean - E_mean * E_mean;

    cmatrix_free(Hpsi);

    double small_ts[] = {1e-4, 1e-3};
    for (size_t k = 0; k < sizeof(small_ts) / sizeof(small_ts[0]); k++) {
      double t = small_ts[k];
      double L_t = ising_loschmidt_echo(eig_f, psi_i, t);
      double L_taylor = 1.0 - t * t * var_Hf;
      char label[128];
      snprintf(label, sizeof label,
               "t=%.0e: L(t) matches short-time Taylor expansion [1 - t^2 * "
               "Var(H_f)]",
               t);
      check_close(L_t, L_taylor, 1e-8, label);
    }

    /* No-quench identity: if \psi_i is (numerically) an eigenvector of H_f
     * itself, L(t)=1 for all t. Use h_f == h_i here so \psi_i is an eigenstate
     * of this H_f, not h_f=2.0 one above. */
    cmatrix_t *Hsame_dense = ising_hamiltonian_dense(N, J, h_i, 1);
    eigen_t *eig_same = cmatrix_eigh_complex(Hsame_dense);
    if (eig_same) {
      double t_tests[] = {0.5, 1.0, 5.0};

      for (size_t k = 0; k < sizeof(t_tests) / sizeof(t_tests[0]); k++) {
        double L_same = ising_loschmidt_echo(eig_same, psi_i, t_tests[k]);
        char label[128];
        snprintf(label, sizeof label,
                 "no-quench (h_f=h_i) t=%.1f: L(t) = 1 for all t", t_tests[k]);

        check_close(L_same, 1.0, 1e-8, label);
      }

      eigen_free(eig_same);
    }

    cmatrix_free(Hsame_dense);
    eigen_free(eig_f);
  }

  cmatrix_free(Hf_dense);
  lanczos_free(gs_i);
  sparse_free(Hi_sparse);
}

int main(void) {
  printf("=== Transverse-Field Ising Model: ED + exact solution tests ===\n");

  test_exact_vs_full_ed();
  test_hermiticity();
  test_z2_parity_reduction();
  test_entanglement_entropy();
  test_quench_dynamics();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_ising_chain checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_ising_chain check(s) FAILED.\n", failures);
    return 1;
  }
}
