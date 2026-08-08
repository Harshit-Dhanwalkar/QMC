/*
 * Multi-Qubit States - Trivial vs. Universal
 *
 * Single-qubit gates alone never produce entanglement (product state, zero
 * entropy), while a single CNOT is enough to build a Bell/GHZ state (maximal
 * entanglement, exponential classical representation cost).
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/qubits.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  printf(" > Multi-Qubit States: Trivial vs. Universal\n\n");

  int n1 = 3;
  cvector_t *psi1 = qstate_alloc(n1);
  for (int q = 0; q < n1; q++) {
    qstate_apply_gate1(psi1, n1, q, hadamard_gate);
  }

  printf("   3 qubits, Hadamard on each (no entangling gate):\n");
  for (int q = 0; q < n1; q++) {
    cmatrix_t *rho = qstate_reduced_density_single(psi1, n1, q);
    double S = von_neumann_entropy_2x2(rho);

    cmatrix_free(rho);

    printf("     qubit %d entanglement entropy = %.6f bits (0 = unentangled)\n",
           q, S);
  }

  cvector_free(psi1);
  printf("\n");

  int n2 = 2;
  cvector_t *psi2 = qstate_alloc(n2);
  qstate_apply_gate1(psi2, n2, 0, hadamard_gate);
  qstate_apply_cnot(psi2, n2, 0, 1);

  printf("   Bell state (H on qubit0, then CNOT(0,1)):\n");
  printf("     P(|00>)=%.4f  P(|01>)=%.4f  P(|10>)=%.4f  P(|11>)=%.4f\n",
         qstate_probability(psi2, 0), qstate_probability(psi2, 1),
         qstate_probability(psi2, 2), qstate_probability(psi2, 3));
  cmatrix_t *rho_bell = qstate_reduced_density_single(psi2, n2, 0);
  double S_bell = von_neumann_entropy_2x2(rho_bell);

  cmatrix_free(rho_bell);

  printf("     entanglement entropy = %.6f bits (1 = maximally entangled)\n\n",
         S_bell);

  cvector_free(psi2);

  int n3 = 3;
  cvector_t *psi3 = qstate_alloc(n3);
  qstate_apply_gate1(psi3, n3, 0, hadamard_gate);
  qstate_apply_cnot(psi3, n3, 0, 1);
  qstate_apply_cnot(psi3, n3, 0, 2);

  printf("   GHZ state (H + CNOT(0,1) + CNOT(0,2), n=3):\n");
  printf("     P(|000>)=%.4f  P(|111>)=%.4f\n", qstate_probability(psi3, 0),
         qstate_probability(psi3, 7));

  cvector_free(psi3);
  printf("\n");

  printf("   Classical state-vector cost (exponential by construction):\n");
  printf("   n_qubits  dimension  memory (complex_t, MB)\n");
  for (int n = 2; n <= 20; n += 4) {
    long long dim = 1LL << n;
    double mb = (double)(dim * sizeof(complex_t)) / (1024.0 * 1024.0);

    printf("   %8d  %9lld  %10.4f\n", n, dim, mb);
  }

  return 0;
}
