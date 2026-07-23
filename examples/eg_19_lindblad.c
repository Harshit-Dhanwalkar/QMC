/*
 * Lindblad Equation - Open Quantum Systems
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/lindblad.h"
#include "../physics/rabi.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  printf(" > Lindblad Master Equation: Open Quantum Systems\n\n");

  // 1. Amplitude damping (T1 relaxation) of an excited qubit
  printf("   Amplitude damping: qubit starts excited (|1><1|), \\gamma=1.0\n");
  printf("   %6s  %10s  %10s  %10s\n", "t", "\\rho11", "\\exp(-t)", "purity");
  {
    double gamma = 1.0, dt = 1e-3;
    cmatrix_t *H = cmatrix_alloc(2, 2);
    cmatrix_t *rho = cmatrix_alloc(2, 2);
    CMAT(rho, 1, 1) = c_real(1.0);
    cmatrix_t *L = lindblad_amplitude_damping_op(1, 0, gamma);
    cmatrix_t *ops[1] = {L};

    for (int checkpoint = 0; checkpoint <= 5; checkpoint++) {
      double t = checkpoint * 0.5;
      if (checkpoint > 0) {
        lindblad_evolve(rho, H, ops, 1, dt, 500); // advance by 0.5 each time
      }
      printf("   %6.2f  %10.6f  %10.6f  %10.6f\n", t, CMAT(rho, 1, 1).re,
             exp(-gamma * t), density_purity(rho));
    }

    cmatrix_free(H);
    cmatrix_free(rho);
    cmatrix_free(L);
  }
  printf("\n");

  // 2. Pure dephasing (T2) of a |+> superposition
  printf("   Pure dephasing: qubit starts in |+> = (|0>+|1>)/\\sqrt2, "
         "\\gamma=1.0\n");
  printf("   %6s  %10s  %10s  %10s\n", "t", "|\\rho01|", "0.5*\\exp(-t)",
         "purity");
  {
    double gamma = 1.0, dt = 1e-3;
    cmatrix_t *H = cmatrix_alloc(2, 2);
    cmatrix_t *rho = cmatrix_alloc(2, 2);
    CMAT(rho, 0, 0) = c_real(0.5);
    CMAT(rho, 0, 1) = c_real(0.5);
    CMAT(rho, 1, 0) = c_real(0.5);
    CMAT(rho, 1, 1) = c_real(0.5);
    cmatrix_t *L = lindblad_dephasing_op(1, 0, gamma);
    cmatrix_t *ops[1] = {L};

    for (int checkpoint = 0; checkpoint <= 5; checkpoint++) {
      double t = checkpoint * 0.5;
      if (checkpoint > 0) {
        lindblad_evolve(rho, H, ops, 1, dt, 500);
      }
      printf("   %6.2f  %10.6f  %10.6f  %10.6f\n", t, c_abs(CMAT(rho, 0, 1)),
             0.5 * exp(-gamma * t), density_purity(rho));
    }

    cmatrix_free(H);
    cmatrix_free(rho);
    cmatrix_free(L);
  }
  printf("\n");

  // 3. Closed-system sanity check: no dissipation reproduces Rabi
  printf("   Closed system (no jump operators) reproduces exact Rabi "
         "oscillation:\n");
  {
    double Omega = 2.0, Delta = 1.0, dt = 1e-4;
    cmatrix_t *H = cmatrix_alloc(2, 2);
    CMAT(H, 0, 0) = c_real(0.5 * Delta);
    CMAT(H, 0, 1) = c_real(0.5 * Omega);
    CMAT(H, 1, 0) = c_real(0.5 * Omega);
    CMAT(H, 1, 1) = c_real(-0.5 * Delta);
    cmatrix_t *rho = cmatrix_alloc(2, 2);
    CMAT(rho, 0, 0) = c_real(1.0);

    printf("   %6s  %14s  %14s\n", "t", "\\rho11 (density)", "P_e (exact)");
    for (int checkpoint = 0; checkpoint <= 4; checkpoint++) {
      double t = checkpoint * 0.5;
      if (checkpoint > 0) {
        lindblad_evolve(rho, H, NULL, 0, dt, 5000);
      }
      printf("   %6.2f  %14.6f  %14.6f\n", t, CMAT(rho, 1, 1).re,
             rabi_excited_probability(t, Omega, Delta));
    }

    cmatrix_free(H);
    cmatrix_free(rho);
  }
  printf("\n");

  // 4. Measurement collapse
  printf(
      "   Measurement collapse of a mixed state (\\rho00=0.3, \\rho11=0.7):\n");
  {
    cmatrix_t *rho = cmatrix_alloc(2, 2);
    CMAT(rho, 0, 0) = c_real(0.3);
    CMAT(rho, 1, 1) = c_real(0.7);
    double u = 0.42; // caller-supplied uniform random draw
    int outcome = density_measure_computational_basis(rho, u);
    printf("     u=%.2f -> collapsed to |%d><%d|\n", u, outcome, outcome);

    cmatrix_free(rho);
  }

  return 0;
}
