/*
 * Variational Quantum Eigensolver
 *
 * Hybrid quantum-classical ground-state solver: a hardware-efficient
 * ansatz circuit (qubits.c) prepares a trial state, and the same
 * golden-section optimizer VMC uses (variational.c) tunes its parameters
 * to minimize the energy expectation value -- the discrete-circuit analog
 * of VMC's continuous trial wavefunction.
 *
 * Demonstrates two cases: an exactly-solvable single qubit (closed-form
 * ground energy), then the transverse-field Ising model on a small open
 * chain, cross-checked against exact diagonalization.
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
  printf(" > Variational Quantum Eigensolver\n\n");

  printf(
      "   Ansatz: n_layers of [RY(\\theta) on every qubit, then CNOT chain]\n");
  printf("   Optimizer: coordinate descent (golden-section line search per "
         "parameter,\n");
  printf("   Jastrow parameter b).\n\n");

  printf("   --- Single qubit: H = a*X + b*Z ---\n");
  double a = 1.3, b = 0.7;
  cmatrix_t *H1 = cmatrix_alloc(2, 2);
  CMAT(H1, 0, 0) = c_real(b);
  CMAT(H1, 0, 1) = c_real(a);
  CMAT(H1, 1, 0) = c_real(a);
  CMAT(H1, 1, 1) = c_real(-b);

  double E1_exact = -sqrt(a * a + b * b);

  double best1 = 1e9;
  for (int trial = 0; trial < 5; trial++) {
    vqe_result_t r = vqe_run(1, 2, H1, 6, M_PI, 1000ULL + (uint64_t)trial);
    if (r.energy < best1) {
      best1 = r.energy;
    }

    free(r.theta_opt);
  }

  printf("   a=%.2f b=%.2f: VQE E=%.6f   exact=-\\sqrt(a^2+b^2)=%.6f\n\n", a, b,
         best1, E1_exact);
  cmatrix_free(H1);

  printf("   --- Transverse-field Ising model, n=3 open chain ---\n");
  printf("   H = -J*\\sum(Z_i Z_i+1) - h*\\sum(X_i)\n\n");
  int n = 3;
  double J = 1.0, h = 0.5;
  cmatrix_t *H2 = vqe_build_tfim(n, J, h);

  cmatrix_t *H2_copy = cmatrix_copy(H2);
  eigen_t *eig = cmatrix_eigh_complex(H2_copy);
  cmatrix_free(H2_copy);
  double E2_exact = eig->eigenvalues[0];
  eigen_free(eig);

  double best2 = 1e9;
  for (int trial = 0; trial < 8; trial++) {
    vqe_result_t r = vqe_run(n, 3, H2, 8, M_PI, 2000ULL + (uint64_t)trial);
    printf("   trial %d: E=%.6f\n", trial, r.energy);
    if (r.energy < best2) {
      best2 = r.energy;
    }

    free(r.theta_opt);
  }

  printf("\n   J=%.2f h=%.2f: VQE (best of 8 restarts) E=%.6f   exact "
         "(diagonalization)=%.6f\n",
         J, h, best2, E2_exact);
  printf("   gap: %.6f Hartree-equivalent units\n\n", best2 - E2_exact);

  cmatrix_free(H2);

  return 0;
}
