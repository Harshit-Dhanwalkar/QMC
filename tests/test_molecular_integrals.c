/*
Test: General Gaussian-type-orbital molecular integrals (McMurchie-Davidson).

1. Boys function: F_0(0)=1 exact identity, plus a couple of known values
   cross-checked against closed-form :
     F_0(x) = \sqrt(\pi / (4 * x)) erf(\sqrt(x))
2. s-s overlap and s-s nuclear attraction against their closed-form
   Gaussian-product-theorem formulas (independent of general McMurchie-Davidson
   code path being tested).
3. molint_normalize_contraction: a normalized basis function's self-overlap must
   be exactly 1.
4. ERI 8-fold permutational symmetry: (ij|kl) must equal every one of its 7
   index-permutation equivalents.
5. The end-to-end validation: full restricted-HF SCF on H2/STO-3G at R=1.4 bohr
   must reproduce the theoretical Szabo & Ostlund result (-1.1167 Hartree) using
   only this integrals engine + a RHF loop written in this test file (not
   reusing hartree_fock.c, which is a radial-grid solver for atoms, not an
   LCAO-Gaussian-basis solver for molecules).
*/

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
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

static void test_boys_function(void) {
  printf("test_boys_function:\n");

  check_close(boys_function(0, 0.0), 1.0, 1e-14, "F_0(0) = 1 exactly");

  for (int i = 0; i < 3; i++) {
    double x = (double[]){0.5, 2.0, 8.0}[i];
    double closed = sqrt(M_PI / (4.0 * x)) * erf(sqrt(x));
    char label[64];

    snprintf(label, sizeof(label), "F_0(%.1f) matches closed form", x);
    check_close(boys_function(0, x), closed, 1e-10, label);
  }

  double F[7];
  boys_function_array(6, 3.0, F);
  check_true(F[0] > F[1] && F[1] > F[2] && F[2] > F[3],
             "boys_function_array is monotonically decreasing in n (F_n(x) "
             "decreasing in n for x>0)");
}

static void test_ss_overlap_and_nuclear_closed_form(void) {
  printf("test_ss_overlap_and_nuclear_closed_form:\n");

  double a = 1.2, b = 0.8;
  double A[3] = {0, 0, 0}, B[3] = {0, 0, 1.4};
  const double alphas_a[1] = {a};
  const double alphas_b[1] = {b};
  const double coeff[1] = {1.0};
  basis_function_t *fa = basis_function_alloc(0, 0, 0, A, 1, alphas_a, coeff);
  basis_function_t *fb = basis_function_alloc(0, 0, 0, B, 1, alphas_b, coeff);

  double S = gto_overlap(fa, fb);
  double p = a + b, R2 = 1.4 * 1.4;
  double S_closed = pow(M_PI / p, 1.5) * exp(-a * b / p * R2);
  check_close(S, S_closed, 1e-12, "s-s overlap matches closed form");

  const double C[3] = {0, 0, 0.7};
  double V = gto_nuclear_attraction(fa, fb, C);
  const double P[3] = {(a * A[0] + b * B[0]) / p, (a * A[1] + b * B[1]) / p,
                       (a * A[2] + b * B[2]) / p};
  double RPC2 = (P[0] - C[0]) * (P[0] - C[0]) + (P[1] - C[1]) * (P[1] - C[1]) +
                (P[2] - C[2]) * (P[2] - C[2]);
  double V_closed =
      2.0 * M_PI / p * exp(-a * b / p * R2) * boys_function(0, p * RPC2);
  check_close(V, V_closed, 1e-12, "s-s nuclear attraction matches closed form");

  basis_function_free(fa);
  basis_function_free(fb);
}

static void test_normalization(void) {
  printf("test_normalization:\n");

  const double center[3] = {0, 0, 0};
  basis_function_t *bf = molint_basis_sto3g_h(center);
  check_true(bf != NULL, "molint_basis_sto3g_h allocates");
  double self_overlap = gto_overlap(bf, bf);
  check_close(self_overlap, 1.0, 1e-12,
              "normalized STO-3G H 1s has self-overlap = 1");
  basis_function_free(bf);
}

static void test_eri_permutational_symmetry(void) {
  printf("test_eri_permutational_symmetry:\n");

  const double A[3] = {0, 0, 0}, B[3] = {0, 0, 1.4};
  basis_function_t *fa = molint_basis_sto3g_h(A);
  basis_function_t *fb = molint_basis_sto3g_h(B);

  double v_abab = gto_eri(fa, fb, fa, fb);
  double v_baab = gto_eri(fb, fa, fa, fb);
  double v_abba = gto_eri(fa, fb, fb, fa);
  double v_abba_swapped = gto_eri(fb, fa, fb, fa);

  check_close(v_abab, v_baab, 1e-12, "(ab|ab) = (ba|ab)");
  check_close(v_abab, v_abba, 1e-12, "(ab|ab) = (ab|ba)");
  check_close(v_abab, v_abba_swapped, 1e-12, "(ab|ab) = (ba|ba)");

  basis_function_free(fa);
  basis_function_free(fb);
}

static void test_molecular_eri_tensor_symmetry(void) {
  printf("test_molecular_eri_tensor_symmetry:\n");

  const double A[3] = {0, 0, 0}, B[3] = {0, 0, 1.4};
  basis_function_t *funcs[2] = {molint_basis_sto3g_h(A),
                                molint_basis_sto3g_h(B)};
  double *eri = molecular_eri_tensor(funcs, 2);
  check_true(eri != NULL, "molecular_eri_tensor allocates");

  check_close(MOLINT_ERI(eri, 2, 0, 1, 0, 1), MOLINT_ERI(eri, 2, 1, 0, 1, 0),
              1e-14, "tensor: (01|01) = (10|10)");
  check_close(MOLINT_ERI(eri, 2, 0, 0, 1, 1), MOLINT_ERI(eri, 2, 1, 1, 0, 0),
              1e-14, "tensor: (00|11) = (11|00)");
  check_close(MOLINT_ERI(eri, 2, 0, 1, 1, 0),
              gto_eri(funcs[0], funcs[1], funcs[1], funcs[0]), 1e-14,
              "tensor: (01|10) matches direct gto_eri call");

  free(eri);

  basis_function_free(funcs[0]);
  basis_function_free(funcs[1]);
}

/* Minimal RHF SCF for a 2-basis-function closed-shell (2 electron) molecule */
static double h2_sto3g_rhf_energy(double R_bond) {
  const double A[3] = {0, 0, 0}, B[3] = {0, 0, R_bond};
  basis_function_t *funcs[2] = {molint_basis_sto3g_h(A),
                                molint_basis_sto3g_h(B)};

  const double charge[2] = {1.0, 1.0};
  const double centers[2][3] = {{0, 0, 0}, {0, 0, R_bond}};
  molecule_t *mol = molecule_alloc(2, charge, centers);

  cmatrix_t *S = molecular_overlap_matrix(funcs, 2);
  cmatrix_t *Hcore = molecular_core_hamiltonian(funcs, 2, mol);
  double *eri = molecular_eri_tensor(funcs, 2);

  /* S^{-1/2} via S's own eigendecomposition. */
  eigen_t *eig_S = cmatrix_eigh_complex(S);
  cmatrix_t *Shalf = cmatrix_alloc(2, 2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      complex_t sum = {0, 0};

      for (int k = 0; k < 2; k++) {
        double inv_sqrt = 1.0 / sqrt(eig_S->eigenvalues[k]);
        complex_t cik = CMAT(eig_S->eigenvectors, i, k);
        complex_t cjk = CMAT(eig_S->eigenvectors, j, k);

        sum.re += cik.re * inv_sqrt * cjk.re;
      }

      CMAT(Shalf, i, j) = sum;
    }
  }

  double D[2][2] = {{0, 0}, {0, 0}};
  double E_elec = 0.0, E_elec_old = 0.0;

  for (int it = 0; it < 50; it++) {
    cmatrix_t *F = cmatrix_alloc(2, 2);

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        double v = CMAT(Hcore, i, j).re;

        for (int k = 0; k < 2; k++) {
          for (int m = 0; m < 2; m++) {
            v += D[k][m] * (MOLINT_ERI(eri, 2, i, j, m, k) -
                            0.5 * MOLINT_ERI(eri, 2, i, k, m, j));
          }
        }

        CMAT(F, i, j) = (complex_t){v, 0.0};
      }
    }

    // F' = Shalf^T F Shalf (Shalf is real symmetric)
    cmatrix_t *Fp = cmatrix_alloc(2, 2);
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        double v = 0.0;

        for (int p = 0; p < 2; p++) {
          for (int q = 0; q < 2; q++) {
            v += CMAT(Shalf, p, i).re * CMAT(F, p, q).re * CMAT(Shalf, q, j).re;
          }
        }

        CMAT(Fp, i, j) = (complex_t){v, 0.0};
      }
    }

    eigen_t *eig_F = cmatrix_eigh_complex(Fp);
    // C = Shalf * C'
    double C[2][2];
    for (int i = 0; i < 2; i++) {
      for (int k = 0; k < 2; k++) {
        double v = 0.0;

        for (int p = 0; p < 2; p++) {
          v += CMAT(Shalf, i, p).re * CMAT(eig_F->eigenvectors, p, k).re;
        }

        C[i][k] = v;
      }
    }

    double D_new[2][2];
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        D_new[i][j] = 2.0 * C[i][0] * C[j][0]; // 1 occupied MO, 2 e-
      }
    }

    E_elec = 0.0;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        E_elec += 0.5 * D_new[i][j] * (CMAT(Hcore, i, j).re + CMAT(F, i, j).re);
      }
    }

    D[0][0] = D_new[0][0];
    D[0][1] = D_new[0][1];
    D[1][0] = D_new[1][0];
    D[1][1] = D_new[1][1];

    cmatrix_free(F);
    cmatrix_free(Fp);
    eigen_free(eig_F);

    if (fabs(E_elec - E_elec_old) < 1e-12) {
      break;
    }

    E_elec_old = E_elec;
  }

  double E_total = E_elec + molecule_nuclear_repulsion(mol);

  eigen_free(eig_S);
  cmatrix_free(S);
  cmatrix_free(Hcore);
  cmatrix_free(Shalf);
  free(eri);
  molecule_free(mol);
  basis_function_free(funcs[0]);
  basis_function_free(funcs[1]);

  return E_total;
}

static void test_h2_sto3g_rhf_benchmark(void) {
  printf("test_h2_sto3g_rhf_benchmark:\n");

  double E = h2_sto3g_rhf_energy(1.4);
  printf("  H2/STO-3G RHF @ R=1.4 bohr: E = %.6f Hartree\n", E);

  /* Reference: Szabo & Ostlund, "Modern Quantum Chemistry", the STO-3G H2
   * benchmark at this exact geometry: -1.1167 Hartree. */
  check_close(E, -1.116714, 1e-5,
              "H2/STO-3G RHF matches validation (-1.116714 Hartree)");
  check_close(E, -1.1167, 2e-4,
              "H2/STO-3G RHF matches theoretical Szabo & Ostlund value");
}

int main(void) {
  test_boys_function();
  test_ss_overlap_and_nuclear_closed_form();
  test_normalization();
  test_eri_permutational_symmetry();
  test_molecular_eri_tensor_symmetry();
  test_h2_sto3g_rhf_benchmark();

  if (failures == 0) {
    printf("\nAll test_molecular_integrals checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_molecular_integrals check(s) FAILED.\n", failures);
    return 1;
  }
}
