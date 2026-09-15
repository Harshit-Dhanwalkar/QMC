/*
 * VQE Beyond TFIM: A General Pauli-String Hamiltonian
 *
 * This example demonstrates the spin-1/2 XXZ Heisenberg chain - same model
 * spin_chain.c diagonalizes exactly via symmetry-sector sparse diagonalization
 * (not dense tensor products) - using the spin-1/2 <-> qubit mapping:
 *   S^z_j = Z_j / 2
 *   S^+_j S^-_k + S^-_j S^+_k = (X_j X_k + Y_j Y_k) / 2
 * so H = \sum_j [Jz S^z_j S^z_{j+1} + (Jxy/2)(S^+_j S^-_{j+1} + h.c.)]
 *      = \sum_j [(Jz/4) Z_j Z_{j+1} + (Jxy/4)(X_j X_{j+1} + Y_j Y_{j+1})]
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../physics/vqe.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > VQE with a General Pauli-String Hamiltonian: the XXZ Chain\n\n");

  int N = 4;
  double Jz = 1.3, Jxy = 0.7;
  printf("  N=%d sites, Jz=%.2f, Jxy=%.2f, periodic boundary conditions\n\n", N,
         Jz, Jxy);

  // Build H = sum_j [(Jz/4) Z_jZ_{j+1} + (Jxy/4)(X_jX_{j+1} + Y_jY_{j+1})]
  // as a Pauli sum: 3 terms per bond (ZZ, XX, YY), N bonds for PBC.
  const char *ops = "ZXY";
  const double base_coeff[3] = {Jz / 4.0, Jxy / 4.0, Jxy / 4.0};

  int n_terms = 3 * N;
  char **strings = malloc((size_t)n_terms * sizeof *strings);
  double *coeffs = malloc((size_t)n_terms * sizeof *coeffs);

  int t = 0;
  for (int j = 0; j < N; j++) {
    int jn = (j + 1) % N;

    for (int o = 0; o < 3; o++) {
      strings[t] = malloc((size_t)N + 1);
      for (int q = 0; q < N; q++) {
        strings[t][q] = (q == j || q == jn) ? ops[o] : 'I';
      }
      strings[t][N] = '\0';
      coeffs[t] = base_coeff[o];
      t++;
    }
  }

  cmatrix_t *H = vqe_build_pauli_hamiltonian(N, (const char *const *)strings,
                                             coeffs, n_terms);

  printf("  Pauli-string decomposition (%d terms):\n", n_terms);
  for (int i = 0; i < n_terms; i++) {
    printf("    %+.4f * %s\n", coeffs[i], strings[i]);
  }
  printf("\n");

  cmatrix_t *H_copy = cmatrix_copy(H);
  eigen_t *eig = cmatrix_eigh_complex(H_copy);
  cmatrix_free(H_copy);
  double E_exact = eig->eigenvalues[0];
  eigen_free(eig);

  printf("  Exact ground energy (dense diagonalization of the Pauli sum): "
         "%.8f\n",
         E_exact);

  double best = 1e9;
  int n_restarts = 10;
  for (int trial = 0; trial < n_restarts; trial++) {
    vqe_result_t r = vqe_run(N, 4, H, 10, M_PI, 3000ULL + (uint64_t)trial);
    if (r.energy < best) {
      best = r.energy;
    }

    free(r.theta_opt);
  }

  printf("  VQE (best of %d restarts): E=%.6f\n", n_restarts, best);
  printf("  gap above exact ground energy: %.6f\n", best - E_exact);

  for (int i = 0; i < n_terms; i++) {
    free(strings[i]);
  }
  free(strings);
  free(coeffs);
  cmatrix_free(H);

  return 0;
}
