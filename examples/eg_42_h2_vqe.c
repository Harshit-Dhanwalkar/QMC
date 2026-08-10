/*
 * H2 VQE: from Atomic-Orbital Integrals to a Quantum Circuit
 *
 * This is fdemo tying together three modules built : molecular_integrals.c (AO
 * integrals), second_quant.c (Jordan-Wigner mapping of molecular electronic
 * Hamiltonian), and vqe.c (variational quantum eigensolver).
 *
 * Pipeline:
 *   1. Build STO-3G AO integrals for H2 (molecular_integrals.c)
 *   2. Restricted Hartree-Fock -> MO coefficients C
 *   3. AO -> MO integral transform
 *   4. Jordan-Wigner: build the 16x16 (4-qubit) second-quantized
 *      Hamiltonian on interleaved spin-orbitals (second_quant.c)
 *   5. Exact diagonalization: the FCI ground state, for comparison
 *   6. VQE: hardware-efficient ansatz, coordinate-descent optimizer
 *      (vqe.c), starting from a random 4-qubit state with NO knowledge
 *      of the answer, converging toward the same ground state via
 *      only energy expectation values -- exactly the classical-optimizer-
 *      in-the-loop workflow a real quantum computer would run, just
 *      simulated classically here.
 *
 * NOTE: global minimum of full unrestricted-particle-number 16-dimensional
 * Hamiltonian sits in physical 2-electron sector for this system, so VQE's
 * hardware-efficient ansatz, which does not itself conserve particle number,
 * can be run directly with no penalty term and still find right answer, which
 * this example checks explicitly.
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

/* RHF returning both electronic energy and MO coefficients C, needed for the
 * AO->MO transform. */
static double rhf(cmatrix_t *Hcore, double *eri, double C_out[2][2],
                  cmatrix_t *S) {
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
  cmatrix_free(Shalf);

  return E_elec;
}

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

int main(void) {
  printf(" > Real H2 VQE: Atomic Integrals -> Jordan-Wigner -> VQE\n\n");

  double R_bond = 1.4;
  double A[3] = {0, 0, 0}, B[3] = {0, 0, R_bond};
  basis_function_t *funcs[2] = {molint_basis_sto3g_h(A),
                                molint_basis_sto3g_h(B)};
  double charge[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R_bond}};
  molecule_t *mol = molecule_alloc(2, charge, centers);

  printf("Step 1: STO-3G AO integrals for H2 at R=%.2f bohr\n\n", R_bond);
  cmatrix_t *S = molecular_overlap_matrix(funcs, 2);
  cmatrix_t *Hcore = molecular_core_hamiltonian(funcs, 2, mol);
  double *eri_ao = molecular_eri_tensor(funcs, 2);
  double nuclear_repulsion = molecule_nuclear_repulsion(mol);

  printf("Step 2: Restricted Hartree-Fock\n\n");
  double C[2][2];
  double E_elec_rhf = rhf(Hcore, eri_ao, C, S);
  double E_rhf = E_elec_rhf + nuclear_repulsion;
  printf("  RHF total energy: %.6f Hartree\n\n", E_rhf);

  printf("Step 3: AO -> MO integral transform\n\n");
  double h_mo[4];
  double *eri_mo = malloc(16 * sizeof(double));
  ao_to_mo(Hcore, eri_ao, C, h_mo, eri_mo);

  printf("Step 4: Jordan-Wigner molecular Hamiltonian "
         "(2 spatial orbitals -> 4 spin-orbitals -> 4 qubits, 16x16)\n\n");
  cmatrix_t *H = second_quant_build_molecular_hamiltonian(2, h_mo, eri_mo,
                                                          nuclear_repulsion);

  printf("Step 5: Exact diagonalization (FCI within this basis)\n\n");
  cmatrix_t *H_copy = cmatrix_alloc(16, 16);
  for (int i = 0; i < 16 * 16; i++) {
    H_copy->data[i] = H->data[i];
  }
  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  double E_fci = eig->eigenvalues[0];

  double N_expect = 0.0;
  for (int state = 0; state < 16; state++) {
    complex_t amp = CMAT(eig->eigenvectors, state, 0);
    double p = amp.re * amp.re + amp.im * amp.im;

    N_expect += p * __builtin_popcount((unsigned)state);
  }
  printf("  FCI ground state: %.6f Hartree  (<N electrons> = %.4f)\n\n", E_fci,
         N_expect);

  printf("Step 6: VQE (4 qubits, hardware-efficient ansatz, coordinate-descent "
         "optimizer)\n\n");
  vqe_result_t vqe_res = vqe_run(4, 3, H, 8, 0.6, 20260810ULL);
  printf("  VQE converged energy: %.6f Hartree\n", vqe_res.energy);
  printf("  Error vs. exact FCI : %.6f Hartree (%.4f mHartree)\n\n",
         vqe_res.energy - E_fci, 1000.0 * (vqe_res.energy - E_fci));

  printf("=== Summary ===\n\n");
  printf("  %-25s %-14.6f\n", "RHF (mean-field)", E_rhf);
  printf("  %-25s %-14.6f\n", "FCI (exact, this basis)", E_fci);
  printf("  %-25s %-14.6f\n", "VQE (variational)", vqe_res.energy);
  printf("\n  Correlation energy captured by FCI/VQE vs RHF: %.6f Hartree\n",
         E_rhf - E_fci);
  printf("  \"Chemical accuracy\" threshold is usually quoted as 1.6 mHartree; "
         "VQE landed %.4f mHartree from exact.\n",
         1000.0 * fabs(vqe_res.energy - E_fci));

  free(vqe_res.theta_opt);
  eigen_free(eig);
  cmatrix_free(H_copy);
  cmatrix_free(H);
  cmatrix_free(S);
  cmatrix_free(Hcore);
  free(eri_ao);
  free(eri_mo);
  molecule_free(mol);
  basis_function_free(funcs[0]);
  basis_function_free(funcs[1]);

  return 0;
}
