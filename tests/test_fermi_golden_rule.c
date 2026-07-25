/*
Test + demonstration: Fermi's Golden Rule, validated against
diagonalization-based time evolution

Setup: couple one initial state |0> (energy 0) to M "quasi-continuum"
final states spread uniformly over an energy window W (so the density of
states rho = M/W = 1/dE, dE = level spacing), each with same coupling
V0 and no coupling among the final states themselves (a "star" / Fano-
Friedrichs-model Hamiltonian).

Fermi's Golden Rule predicts P(t) ~ \exp(-\Gamma*t), \Gamma = 2 * \pi * V0^2 *
\rho, valid at intermediate times: t >> 1/W (needed for golden-rule energy
coarse-graining to apply, at very short times decay is quadratic/Zeno-like, not
exponential) and t << 1/dE (before finite-size recurrences from discrete
spectrum set in). The test checks agreement in intermediate window, not at very
short t.
*/

#include "../core/complex.h"
#include "../core/linalg/eigen_generic.h"
#include "../core/matrix.h"
#include "../physics/perturbation.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.6f expected=%.6f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

int main(void) {
  int M = 150;
  double W = 4.0;
  double V0 = 0.05;
  int dim = M + 1;

  double dE = W / (M - 1);
  double rho = 1.0 / dE;

  // Build exact (M+1)x(M+1) Hamiltonian: state 0 = initial (E=0),
  // states 1..M = manifold, each coupled to state 0 with strength V0.
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

  // fermi_golden_rate() itself, called on the same coupling.
  cmatrix_t *Vpert = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    Vpert->data[i] = c_zero();
  }
  CMAT(Vpert, 1, 0) = c_real(V0);
  double Gamma = fermi_golden_rate(Vpert, 0, 1, rho);

  int fail =
      check_close(Gamma, 2.0 * M_PI * V0 * V0 * rho, 1e-10,
                  "fermi_golden_rate() vs hand formula 2*\\pi*V0^2*\\rho");

  // Exact time evolution via diagonalization.
  eigen_t *eig = cmatrix_eigh_generic(H);
  cmatrix_free(H);
  cmatrix_free(Vpert);
  if (!eig) {
    printf("  cmatrix_eigh_generic failed\n");

    return 1;
  }

  double *c = malloc(dim * sizeof(double));
  for (int k = 0; k < dim; k++) {
    c[k] = CMAT(eig->eigenvectors, 0, k).re; // <eigvec_k|0>
  }

  printf("  Gamma (golden rule) = %.6f\n", Gamma);
  printf("  Survival probability vs exp(-Gamma*t), intermediate-time "
         "window:\n");
  double test_times[] = {2.0, 2.5, 3.0, 3.5, 4.0};

  for (int ti = 0; ti < 5; ti++) {
    double t = test_times[ti];
    double re = 0.0, im = 0.0;
    for (int k = 0; k < dim; k++) {
      double phase = -eig->eigenvalues[k] * t;
      re += c[k] * c[k] * cos(phase);
      im += c[k] * c[k] * sin(phase);
    }

    double P = re * re + im * im;
    char label[32];
    snprintf(label, sizeof label, "P_survival(t=%.1f)", t);
    // relative tolerance: Compares an exact many-level simulation against an
    // asymptotic (intermediate-time) formula
    fail |=
        check_close(P, exp(-Gamma * t), 0.05 * exp(-Gamma * t) + 0.02, label);
  }

  free(c);
  eigen_free(eig);

  if (fail) {
    printf("FAILED\n");
    return 1;
  }
  printf("PASS\n");

  return 0;
}
