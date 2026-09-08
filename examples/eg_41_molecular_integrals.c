/*
 * General Molecular Integrals (McMurchie-Davidson) : H2/STO-3G Restricted
 * Hartree-Fock
 *
 * HACK: physics/molecular_integrals.c is a general Gaussian-type-orbital (GTO)
 * integrals engine (arbitrary angular momentum, arbitrary contraction), built
 * specifically so second_quant.c + vqe.c can eventually tackle a real molecule
 * instead of only toy fermionic-lattice Hamiltonians
 * second_quant_build_hopping_hamiltonian() provides.
 *
 * This example runs complete pipeline end to end on the smallest real molecule,
 * H2, in smallest real basis, STO-3G: build the basis, compute
 * overlap/kinetic/nuclear/two-electron integrals, run a restricted-HF SCF, and
 * report total energy.
 * The Theoretical reference result for exactly this system/basis/geometry
 * (Reference: Szabo & Ostlund, "Modern Quantum Chemistry") is -1.1167 Hartree
 * at R=1.4 bohr : this example should land within a few 1e-5 Hartree of that.
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double h2_sto3g_rhf(double R_bond, int verbose) {
  double A[3] = {0, 0, 0}, B[3] = {0, 0, R_bond};
  basis_function_t *funcs[2] = {molint_basis_sto3g_h(A),
                                molint_basis_sto3g_h(B)};

  const double charge[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R_bond}};
  molecule_t *mol = molecule_alloc(2, charge, centers);

  cmatrix_t *S = molecular_overlap_matrix(funcs, 2);
  cmatrix_t *Hcore = molecular_core_hamiltonian(funcs, 2, mol);
  double *eri = molecular_eri_tensor(funcs, 2);

  if (verbose) {
    printf("  Overlap matrix S:\n");
    printf("    [ %.6f  %.6f ]\n", CMAT(S, 0, 0).re, CMAT(S, 0, 1).re);
    printf("    [ %.6f  %.6f ]\n", CMAT(S, 1, 0).re, CMAT(S, 1, 1).re);
    printf("  Core Hamiltonian h = T - \\sum_A Z_A * V_A:\n");
    printf("    [ %.6f  %.6f ]\n", CMAT(Hcore, 0, 0).re, CMAT(Hcore, 0, 1).re);
    printf("    [ %.6f  %.6f ]\n", CMAT(Hcore, 1, 0).re, CMAT(Hcore, 1, 1).re);
    printf("  Unique two-electron integrals:\n");
    printf("    (00|00) = %.6f   (00|11) = %.6f\n",
           MOLINT_ERI(eri, 2, 0, 0, 0, 0), MOLINT_ERI(eri, 2, 0, 0, 1, 1));
    printf("    (01|01) = %.6f   (00|01) = %.6f\n",
           MOLINT_ERI(eri, 2, 0, 1, 0, 1), MOLINT_ERI(eri, 2, 0, 0, 0, 1));
  }

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
  int it;

  for (it = 0; it < 50; it++) {
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
        D_new[i][j] = 2.0 * C[i][0] * C[j][0];
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

  double E_nuc = molecule_nuclear_repulsion(mol);
  double E_total = E_elec + E_nuc;

  if (verbose) {
    printf("  Converged in %d SCF iterations\n", it + 1);
    printf("  E_electronic = %.6f Hartree\n", E_elec);
    printf("  E_nuclear    = %.6f Hartree\n", E_nuc);
  }

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

int main(void) {
  printf(
      " > General Molecular Integrals: H2/STO-3G Restricted Hartree-Fock\n\n");

  printf("  === H2 at R=1.4 bohr (Szabo & Ostlund benchmark geometry) ===\n\n");
  double E = h2_sto3g_rhf(1.4, 1);
  printf(
      "\n  E_total = %.6f Hartree (Theoretical reference: -1.1167 Hartree)\n\n",
      E);

  printf("  === Bond-length scan ===\n\n");
  printf("  %-10s %-12s\n", "R (bohr)", "E (Hartree)");
  for (int i = 0; i < 12; i++) {
    double R = 0.8 + i * 0.15;
    double Escan = h2_sto3g_rhf(R, 0);

    printf("  %-10.3f %-12.6f\n", R, Escan);
  }

  return 0;
}
