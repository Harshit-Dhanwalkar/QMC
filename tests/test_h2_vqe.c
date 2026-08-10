/*
Test: full H2/STO-3G VQE pipeline, end to end :
AO integrals (molecular_integrals.c) -> RHF -> AO-to-MO transform ->
second-quantized Jordan-Wigner Hamiltonian (second_quant.c) -> exact
diagonalization (FCI) -> VQE (vqe.c).
*/

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
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

/* Same RHF as test_molecular_integrals.c, but also returns converged MO
 * coefficient matrix C, needed for AO->MO integral transform below. */
static double h2_sto3g_rhf_with_C(double R_bond, basis_function_t **funcs,
                                  cmatrix_t *Hcore, double *eri,
                                  double C_out[2][2]) {
  cmatrix_t *S = molecular_overlap_matrix(funcs, 2);
  eigen_t *eig_S = cmatrix_eigh_complex(S);
  cmatrix_t *Shalf = cmatrix_alloc(2, 2);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      double sum = 0.0;

      for (int k = 0; k < 2; k++) {
        double inv_sqrt = 1.0 / sqrt(eig_S->eigenvalues[k]);
        sum += CMAT(eig_S->eigenvectors, i, k).re * inv_sqrt *
               CMAT(eig_S->eigenvectors, j, k).re;
      }

      CMAT(Shalf, i, j) = (complex_t){sum, 0.0};
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
    for (int i = 0; i < 2; i++) {
      for (int k = 0; k < 2; k++) {
        double v = 0.0;

        for (int p = 0; p < 2; p++) {
          v += CMAT(Shalf, i, p).re * CMAT(eig_F->eigenvectors, p, k).re;
        }

        C_out[i][k] = v;
      }
    }

    double D_new[2][2];
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        D_new[i][j] = 2.0 * C_out[i][0] * C_out[j][0];
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

    if (fabs(E_elec - E_elec_old) < 1e-13) {
      break;
    }

    E_elec_old = E_elec;
  }

  eigen_free(eig_S);
  cmatrix_free(S);
  cmatrix_free(Shalf);

  (void)R_bond;
  return E_elec;
}

/* AO (2x2 / 2x2x2x2) -> MO transform via the RHF coefficient matrix C, chemist
 * notation preserved throughout. */
static void ao_to_mo(const cmatrix_t *h_ao, const double *eri_ao,
                     double C[2][2], double h_mo[4], double *eri_mo) {
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      double v = 0.0;

      for (int p = 0; p < 2; p++) {
        for (int q = 0; q < 2; q++) {
          v += C[p][i] * CMAT(h_ao, p, q).re * C[q][j];
        }
      }

      h_mo[i * 2 + j] = v;
    }
  }

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        for (int l = 0; l < 2; l++) {
          double v = 0.0;

          for (int p = 0; p < 2; p++) {
            for (int q = 0; q < 2; q++) {
              for (int r = 0; r < 2; r++) {
                for (int s = 0; s < 2; s++) {
                  v += C[p][i] * C[q][j] * C[r][k] * C[s][l] *
                       MOLINT_ERI(eri_ao, 2, p, q, r, s);
                }
              }
            }
          }

          MOLINT_ERI(eri_mo, 2, i, j, k, l) = v;
        }
      }
    }
  }
}

static void test_h2_fci_and_vqe(void) {
  printf("test_h2_fci_and_vqe:\n");

  double R_bond = 1.4;
  double A[3] = {0, 0, 0}, B[3] = {0, 0, R_bond};
  basis_function_t *funcs[2] = {molint_basis_sto3g_h(A),
                                molint_basis_sto3g_h(B)};

  double charge[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R_bond}};
  molecule_t *mol = molecule_alloc(2, charge, centers);

  cmatrix_t *Hcore = molecular_core_hamiltonian(funcs, 2, mol);
  double *eri_ao = molecular_eri_tensor(funcs, 2);
  double nuclear_repulsion = molecule_nuclear_repulsion(mol);

  double C[2][2];
  double E_elec_rhf = h2_sto3g_rhf_with_C(R_bond, funcs, Hcore, eri_ao, C);
  double E_rhf = E_elec_rhf + nuclear_repulsion;
  printf("  RHF total energy: %.6f Hartree\n", E_rhf);
  check_close(E_rhf, -1.116714, 1e-4, "RHF matches earlier validated value");

  double h_mo[4];
  double *eri_mo = malloc(16 * sizeof(double));
  ao_to_mo(Hcore, eri_ao, C, h_mo, eri_mo);

  cmatrix_t *H_molecular = second_quant_build_molecular_hamiltonian(
      2, h_mo, eri_mo, nuclear_repulsion);
  check_true(H_molecular != NULL,
             "second_quant_build_molecular_hamiltonian allocates");
  check_true(H_molecular->nrows == 16,
             "H is 16x16 (2 spatial orbitals -> 4 spin-orbitals -> 2^4)");

  // Exact diagonalization (FCI within this basis)
  cmatrix_t *H_copy = cmatrix_alloc(16, 16);
  for (int i = 0; i < 16 * 16; i++) {
    H_copy->data[i] = H_molecular->data[i];
  }
  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  double E_fci = eig->eigenvalues[0];
  printf("  FCI (exact diagonalization) ground state: %.6f Hartree\n", E_fci);
  check_close(E_fci, -1.137276, 1e-4,
              "FCI matches Python cross-validation (-1.137276 Hartree, also "
              "the well-known literature FCI/STO-3G H2 value)");
  check_true(E_fci < E_rhf - 1e-6,
             "FCI is below RHF (captures correlation energy, as it must "
             "variationally)");

  // NOTE: Ground-state particle number: sum over computational-basis states of
  // |amplitude|^2 * popcount(state), should be exactly 2 (global
  // unrestricted-Fock-space minimum happens to sit in physical 2-electron
  // sector for this system)
  double N_expect = 0.0;
  for (int state = 0; state < 16; state++) {
    complex_t amp = CMAT(eig->eigenvectors, state, 0);
    double p = amp.re * amp.re + amp.im * amp.im;
    int nbits = __builtin_popcount((unsigned)state);

    N_expect += p * nbits;
  }
  check_close(N_expect, 2.0, 1e-6,
              "ground state electron number expectation = 2 exactly");

  eigen_free(eig);
  cmatrix_free(H_copy);

  // NOTE: VQE: unconstrained hardware-efficient ansatz on same full Hamiltonian
  // (no particle-number restriction/penalty). Should converge close to the FCI
  // value.
  vqe_result_t vqe_res = vqe_run(4, 3, H_molecular, 6, 0.6, 20260810ULL);
  printf("  VQE (4 qubits, 3 layers): E = %.6f Hartree (FCI = %.6f)\n",
         vqe_res.energy, E_fci);
  check_true(vqe_res.theta_opt != NULL, "VQE ran (theta_opt allocated)");
  check_close(vqe_res.energy, E_fci, 0.05,
              "VQE converges close to FCI ground state");
  check_true(vqe_res.energy >= E_fci - 1e-6,
             "VQE energy is variational (>= exact ground state)");

  free(vqe_res.theta_opt);
  cmatrix_free(H_molecular);
  cmatrix_free(Hcore);
  free(eri_ao);
  free(eri_mo);
  molecule_free(mol);
  basis_function_free(funcs[0]);
  basis_function_free(funcs[1]);
}

int main(void) {
  test_h2_fci_and_vqe();

  if (failures == 0) {
    printf("\nAll test_h2_vqe checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_h2_vqe check(s) FAILED.\n", failures);
    return 1;
  }
}
