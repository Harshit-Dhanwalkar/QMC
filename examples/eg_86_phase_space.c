/*
 * Kerr cat states in phase space
 *
 * NOTE: A coherent state |alpha = 2.5> evolves under Kerr Hamiltonian
 * H = \chi n (n - 1). Its Wigner function smears into a ring (photon number is
 * conserved, phase is not), and at \chi(t) = \pi / 2 it is exactly Schroedinger
 * cat (\exp(-i \pi/4) |i \alpha> + \exp(i \pi/4) |-i \alpha>) / \sqrt(2): two
 * blobs plus interference fringes that go negative. Negative volume of W
 * measures how non-classical state is
 */

#include "../physics/phase_space.h"
#include "complex.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NMAX 100
#define NX 128

int main(void) {
  const double alpha = 2.5;
  const double dx = 0.125;
  complex_t c[NMAX + 1];
  complex_t psi[NX];
  double xs[NX], ps[NX];
  static double W[NX * NX];

  for (int i = 0; i < NX; i++) {
    xs[i] = ps[i] = (i - (NX - 1) / 2.0) * dx;
  }

  printf(
      " > Kerr evolution of a coherent state, \\alpha = %.1f, <n> = %.2f\n",
      alpha, alpha * alpha);
  printf("  %8s  %10s  %14s  %12s\n", "\\chi(t)", "<n>", "W minimum",
         "neg. volume");

  for (int i = 0; i <= 8; i++) {
    double chit = M_PI / 4.0 * i;

    ps_coherent_state(c, NMAX, alpha, 0.0);
    ps_kerr_evolve(c, NMAX, chit);
    ps_wavefunction(c, NMAX, xs, NX, psi);
    ps_wigner(psi, NX, dx, ps, NX, W);

    double wmin = 0.0;

    for (int k = 0; k < NX * NX; k++) {
      wmin = fmin(wmin, W[k]);
    }

    printf("  %8.4f  %10.5f  %14.6f  %12.6f\n", chit,
           ps_mean_photon_number(c, NMAX), wmin,
           ps_wigner_negativity(W, NX, NX, dx, dx));
  }

  return 0;
}
