/*
 * Verifies potential functions in physics/potentials.c compute
 * correct values at specific points against known analytic results.
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

  // Infinite square well: V=0 for |x|<a, V=\inf elsewhere
  printf("   Infinite square well (a=1.0):\n");
  double a_iw = 1.0;
  potential_array(x, N, V_infinite_well, &a_iw, V);
  int i_mid = N / 2;                     // x=0 (inside)
  int i_wall = (int)((1.5 - (-L)) / dx); // x=1.5 (outside)
  pass &= check("V_infinite_well(0)", V[i_mid], 0.0, 1e-12);
  if (V[i_wall] <= 1e8) {
    printf("   FAIL: V(1.5) should be large (got %.2e)\n", V[i_wall]);
    pass = 0;
  } else {
    printf("   V(0)=0, V(1.5)=inf (outside well): PASS\n");
  }

  // Harmonic oscillator: V(x) = 0.5*\omega^2*x^2
  printf("   Harmonic oscillator (omega=1):\n");
  double omega = 1.0;
  potential_array(x, N, V_harmonic, &omega, V);
  int idx1 = (int)((1.0 - (-L)) / dx);
  int idx2 = (int)((2.0 - (-L)) / dx);
  pass &= check("V_HO(1.0)", V[idx1], 0.5, 0.01);
  pass &= check("V_HO(2.0)", V[idx2], 2.0, 0.01);
  printf("   V(1.0)=%.4f (exp 0.5), V(2.0)=%.4f (exp 2.0)\n", V[idx1], V[idx2]);

  // Step potential: V=0 for x<0, V=V0 for x>=0
  printf("   Step potential (V0=5):\n");
  double V0_step = 5.0;
  potential_array(x, N, V_step, &V0_step, V);
  int mid = N / 2;
  pass &= check("V_step(x<0)", V[mid - 20], 0.0, 1e-12);
  pass &= check("V_step(x>0)", V[mid + 20], 5.0, 1e-12);
  printf("   V(x<0)=%.1f, V(x>0)=%.1f\n", V[mid - 20], V[mid + 20]);

  // Rectangular barrier: V=V0 for 0<x<a, V=0 elsewhere
  printf("   Rectangular barrier (a=1.0, V0=3.0):\n");
  double barrier_params[2] = {1.0, 3.0}; // {a, V0}
  potential_array(x, N, V_barrier, barrier_params, V);
  int i_in = (int)((0.5 - (-L)) / dx);  // inside (0, a)
  int i_out = (int)((1.5 - (-L)) / dx); // outside
  pass &= check("V_barrier(0.5)", V[i_in], 3.0, 1e-12);
  pass &= check("V_barrier(1.5)", V[i_out], 0.0, 1e-12);
  printf("   V(0.5)=%.1f (exp 3.0), V(1.5)=%.1f (exp 0.0)\n", V[i_in],
         V[i_out]);

  // Analytic energy levels
  printf("   Analytic energy levels (hbar=m=1):\n");
  // Infinite well of half-width a: E_n = n^2 \pi^2 / (8 a^2)
  double a_energy = 1.0;
  double E1_well = (1.0 * M_PI * M_PI) / (8.0 * a_energy * a_energy);
  double E1_expected = M_PI * M_PI / 8.0;
  pass &= check("E1_infinite_well", E1_well, E1_expected, 1e-10);
  // Harmonic oscillator: E_n = (n + 1/2) \hbar \omega
  double E0_ho = (0 + 0.5) * omega;
  double E1_ho = (1 + 0.5) * omega;
  pass &= check("E0_harmonic", E0_ho, 0.5, 1e-10);
  pass &= check("E1_harmonic", E1_ho, 1.5, 1e-10);
  printf("   E1_well=%.4f (exp %.4f), E0_HO=%.4f, E1_HO=%.4f\n", E1_well,
         E1_expected, E0_ho, E1_ho);

  free(x);
  free(V);
  if (pass)
    printf("   Potentials test passed.\n");
  else
    printf("   Potentials test FAILED.\n");
  return pass ? 0 : 1;
}
