/*
 * Quantum kicked rotor: dynamical localisation
 *
 * NOTE: The classical standard map at K = 5 is chaotic and the momentum
 * diffuses, <p^2> = D n with D ~ K^2/2. The quantum rotor with hbar = 1 follows
 * it for a few kicks and then stops: <m^2> saturates (dynamical localisation).
 * Prints both in units of hbar^2, then the exact quantum resonance hbar = 4 pi,
 * where <m^2> = (n K / hbar)^2 / 2 grows ballistically.
 */

#include "../physics/kicked_rotor.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { N = 4096, NP = 20000 };

static double g_theta[NP];
static double g_p[NP];

int main(void) {
  const double k = 5.0;
  const double hbar = 1.0;
  cvector_t *psi = cvector_alloc(N);

  if (!psi) {
    return 1;
  }
  krotor_init(psi);
  krotor_classical_init(g_theta, g_p, NP, hbar, 7);
  printf(" > Kicked rotor, K = %.1f, hbar = %.1f\n\n", k, hbar);
  printf("  %6s  %14s  %14s\n", "kicks", "quantum <m^2>", "classical <m^2>");
  int t = 0;

  for (int target = 5; target <= 640; target *= 2) {
    krotor_step(psi, k, hbar, target - t);
    krotor_classical_step(g_theta, g_p, NP, k, target - t);
    t = target;
    printf("  %6d  %14.2f  %14.2f\n", t, krotor_momentum2(psi),
           krotor_classical_p2(g_p, NP) / (hbar * hbar));
  }
  printf("\n  Quantum resonance, hbar = 4 pi:\n");
  krotor_init(psi);
  for (int n = 1; n <= 5; n++) {
    krotor_step(psi, k, 4.0 * M_PI, 1);
    double x = n * k / (4.0 * M_PI);

    printf("  n = %d  <m^2> = %.6f  exact %.6f\n", n, krotor_momentum2(psi),
           0.5 * x * x);
  }
  cvector_free(psi);

  return 0;
}
