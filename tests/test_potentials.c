/*
 * Verifies potential functions in physics/potentials.c compute
 * correct values at specific points with known analytic answers.
 */

#include "../core/utils.h"
#include "../physics/potentials.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check(const char *name, double got, double expected, double tol) {
  int pass = fabs(got - expected) < tol;
  if (!pass)
    printf("   FAIL %s: got %.6e expected %.6e (tol %.1e)\n", name, got,
           expected, tol);
  return pass;
}

int main(void) {
  printf(" > Testing potential functions...\n");
  int pass = 1;

  int N = 201;
  double L = 2.0;
  double dx = 2.0 * L / (N - 1);
  double *x = linspace(-L, L, N);
  double *V = malloc(N * sizeof *V);
  if (!x || !V) {
    printf("FAIL: memory\n");
    return 1;
  }

  // Infinite square well
  printf("   Infinite square well (L=1, centred at 0.5):\n");
  double *x_well = linspace(0.0, 1.0, N);
  potential_infinite_well(V, x_well, N, 1.0);
  // Inside: V=0
  pass &= check("V(0.5)", V[N / 2], 0.0, 1e-12);
  // At wall (x=0): V should be large
  pass &= (V[0] > 1e8);
  if (V[0] <= 1e8)
    printf("   FAIL: V(0) should be large (got %.2e)\n", V[0]);
  else
    printf("   V(0.5)=0, V(0)=inf (walls): PASS\n");
  free(x_well);

  // Harmonic oscillator
  printf("   Harmonic oscillator (m=1, \\omega=1):\n");
  potential_harmonic(V, x, N, 1.0, 1.0);
  // V(x) = ½x^2, check at x=1: V=0.5, x=2: V=2.0
  // Find index closest to x=1
  int idx1 = (int)((1.0 - (-L)) / dx);
  int idx2 = (int)((2.0 - (-L)) / dx);
  pass &= check("V_HO(1.0)", V[idx1], 0.5, 0.01);
  pass &= check("V_HO(2.0)", V[idx2], 2.0, 0.01);
  printf("   V(1.0)=%.4f (exp 0.5), V(2.0)=%.4f (exp 2.0): %s\n", V[idx1],
         V[idx2],
         (fabs(V[idx1] - 0.5) < 0.01 && fabs(V[idx2] - 2.0) < 0.01) ? "PASS"
                                                                    : "FAIL");

  // Step potential
  printf("   Step potential (V0=5, step at x=0):\n");
  potential_step(V, x, N, 5.0);
  int mid = N / 2;
  // Left of step: V=0, right: V=5
  pass &= check("V_step(x<0)", V[mid - 20], 0.0, 1e-12);
  pass &= check("V_step(x>0)", V[mid + 20], 5.0, 1e-12);
  printf("   V(x<0)=%.1f, V(x>0)=%.1f: %s\n", V[mid - 20], V[mid + 20],
         (fabs(V[mid - 20]) < 1e-10 && fabs(V[mid + 20] - 5.0) < 1e-10)
             ? "PASS"
             : "FAIL");

  // Rectangular barrier
  printf("   Rectangular barrier (V0=3, x∈[-0.5,0.5]):\n");
  potential_barrier(V, x, N, -0.5, 0.5, 3.0);
  int i_in = (int)((0.0 - (-L)) / dx);
  int i_out = (int)((1.0 - (-L)) / dx);
  pass &= check("V_barrier(0)", V[i_in], 3.0, 1e-12);
  pass &= check("V_barrier(1.0)", V[i_out], 0.0, 1e-12);
  printf("   V(0)=%.1f (exp 3.0), V(1.0)=%.1f (exp 0.0): %s\n", V[i_in],
         V[i_out],
         (fabs(V[i_in] - 3.0) < 1e-10 && fabs(V[i_out]) < 1e-10) ? "PASS"
                                                                 : "FAIL");

  // Analytic energy levels
  printf("   Analytic energy functions:\n");
  // Infinite well: E_n = n^2\pi ^2/(2L^2), L=1, \hbar=m=1
  double E1_well = energy_infinite_well(1, 1.0, 0.5, 1.0);
  double E1_expected = M_PI * M_PI / 2.0;
  pass &= check("E1_infinite_well", E1_well, E1_expected, 1e-10);
  // HO: E_n = (n+0.5) \hbar \omega
  double E0_ho = energy_harmonic(0, 1.0, 1.0);
  pass &= check("E0_harmonic", E0_ho, 0.5, 1e-10);
  double E1_ho = energy_harmonic(1, 1.0, 1.0);
  pass &= check("E1_harmonic", E1_ho, 1.5, 1e-10);
  printf("   E1_well=%.4f (exp %.4f), E0_HO=%.4f, E1_HO=%.4f: %s\n", E1_well,
         E1_expected, E0_ho, E1_ho, pass ? "PASS" : "FAIL");

  free(x);
  free(V);
  if (pass)
    printf("   Potentials test passed.\n");
  return pass ? 0 : 1;
}
