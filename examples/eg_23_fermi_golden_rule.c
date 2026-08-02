/*
 * Fermi's Golden Rule - From Formula to Exponential Decay
 *
 * fermi_golden_rate() in perturbation.c computes a single number, \Gamma = 2 *
 * \pi * |V_fi|^2 * \rho(E) but point of Golden Rule is that this rate governs
 * genuine exponential decay out of initial state when it's coupled to dense
 * manifold of final states. This example builds that manifold explicitly (a
 * "star" Hamiltonian: one initial state coupled equally to M near-continuum
 * final states, no coupling among final states), evolves it exactly via
 * diagonalization, and shows survival probability settling onto \exp(-\Gamma *
 * t) at intermediate times, after initial quadratic ("quantum Zeno") transient
 * and before finite-size recurrences.
 */

#include "../core/complex.h"
#include "../core/linalg/eigen_generic.h"
#include "../core/matrix.h"
#include "../physics/perturbation.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Fermi's Golden Rule: Coupling to a Quasi-Continuum\n\n");

  int M = 150;      // no. of final states in the manifold
  double W = 4.0;   // energy window the manifold spans
  double V0 = 0.05; // coupling strength, initial state to each final state
  int dim = M + 1;

  double dE = W / (M - 1);
  double rho = 1.0 / dE; // density of states

  cmatrix_t *H = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    H->data[i] = c_zero();
  }

  for (int k = 1; k <= M; k++) {
    double Ek = -W / 2.0 + (k - 1) * dE;
    CMAT(H, k, k) = c_real(Ek);
    CMAT(H, 0, k) = c_real(V0);
    CMAT(H, k, 0) = c_real(V0);
  }

  cmatrix_t *Vpert = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    Vpert->data[i] = c_zero();
  }

  CMAT(Vpert, 1, 0) = c_real(V0);
  double Gamma = fermi_golden_rate(Vpert, 0, 1, rho);
  cmatrix_free(Vpert);

  printf("   Manifold: %d final states over energy window W=%.1f (\\rho=%.2f, "
         "V_0=%.3f)\n",
         M, W, rho, V0);
  printf("   Fermi's Golden Rule: \\Gamma = 2*\\pi * V_0^2 * \\rho = %.6f\n\n",
         Gamma);

  eigen_t *eig = cmatrix_eigh_generic(H);
  cmatrix_free(H);

  double *c = malloc(dim * sizeof(double));
  for (int k = 0; k < dim; k++) {
    c[k] = CMAT(eig->eigenvectors, 0, k).re;
  }

  printf("   Survival probability P(t) = |<initial|\\psi(t)>|^2:\n");
  printf("   %6s  %14s  %14s\n", "t", "exact P(t)", "\\exp(-\\Gamma*t)");
  for (double t = 0.25; t <= 5.0; t += 0.25) {
    double re = 0.0, im = 0.0;
    for (int k = 0; k < dim; k++) {
      double phase = -eig->eigenvalues[k] * t;
      re += c[k] * c[k] * cos(phase);
      im += c[k] * c[k] * sin(phase);
    }
    double P = re * re + im * im;
    printf("   %6.2f  %14.6f  %14.6f\n", t, P, exp(-Gamma * t));
  }
  printf("\n  (short t: quadratic/Zeno transient, not exponential : expected; "
         "t ~ 2-5: settles onto golden-rule exponential decay)\n");

  free(c);
  eigen_free(eig);

  return 0;
}
