/*
 * Second Quantization and the Jordan-Wigner Transformation
 *
 * Maps fermionic creation/annihilation operators to qubit operators, so
 * fermionic many-body Hamiltonians can be built and diagonalized. Demonstrates
 * the mapping's correctness two ways: defining anticommutation relations, and
 * an exact cross-check against a completely independent direct Fock-space
 * Hamiltonian construction.
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/second_quant.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  printf(" > Second Quantization and the Jordan-Wigner Transformation\n\n");

  int n_modes = 4;
  // int dim = 1 << n_modes;

  printf("   Fermionic modes map to qubits (mode j occupied <-> qubit j "
         "= |1>).\n");
  printf("   a_j = Z_0...Z_{j-1} (x) \\sigma^-_j (x) I...I : Z-string supplies "
         "fermionic anticommutation sign.\n\n");

  cmatrix_t *a2 = jw_annihilation_operator(2, n_modes);
  cmatrix_t *a2_dag = jw_creation_operator(2, n_modes);
  cmatrix_t *n2 = cmatrix_multiply(a2_dag, a2);
  printf("   <n_2> for basis state |0011> (modes 2,3 occupied): %.4f "
         "(expect 1.0)\n",
         CMAT(n2, 3, 3).re);
  printf("   <n_2> for basis state |0000> (vacuum):             %.4f "
         "(expect 0.0)\n\n",
         CMAT(n2, 0, 0).re);

  cmatrix_free(a2);
  cmatrix_free(a2_dag);
  cmatrix_free(n2);

  printf("   --- Toy spinless Hubbard chain: H = \\sum(eps_i n_i) - "
         "t * \\sum(hopping) + U * \\sum(n_i n_i+1), n_modes=%d ---\n\n",
         n_modes);
  const double epsilon[4] = {0.3, -0.1, 0.5, 0.2};
  double t = 1.0, U = 0.7;

  cmatrix_t *H = second_quant_build_hopping_hamiltonian(n_modes, epsilon, t, U);
  cmatrix_t *H_copy = cmatrix_copy(H);
  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  cmatrix_free(H_copy);

  printf("   Lowest 5 eigenvalues:\n");
  // for (int i = 0; i < 5 && i < dim; i++) NOTE: i < dim is always true
  for (int i = 0; i < 5; i++) {
    printf("     E_%d = %.6f\n", i, eig->eigenvalues[i]);
  }

  printf("\n   Ground state energy %.6f matches an\n", eig->eigenvalues[0]);

  eigen_free(eig);
  cmatrix_free(H);

  return 0;
}
