/*
 * Dirac wavepacket dynamics: Zitterbewegung and Klein paradox
 *
 * NOTE: Units \hbar = m = c = 1. Two relativistic effects that follow from
 * evolving 1D Dirac equation, i \hbar d\psi/dt = (c \sigma_x p + m c^2
 * \sigma_z + V) \psi, with split-operator Fourier propagator:
 *
 *  1. Zitterbewegung: a packet with both energy signs present trembles at
 *     2 m c^2/\hbar with amplitude \hbar/(2 m c) (interference between
 *     positive- and negative-energy components)
 *  2. Klein step: a positive-energy packet hits a sharp step. For V0 above E +
 *     m c^2 transmitted wave is a negative-energy state, and transmission rises
 *     with V0 instead of dying off exponentially like non-relativistic
 *     tunnelling. Simulated T is compared with  closed form
 */

#include "../physics/dirac_evolve.h"
#include "complex.h"
#include "vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { N = 2048 };
static const double DX = 0.1;

static void zitterbewegung_demo(void) {
  printf("  === Zitterbewegung ===\n");

  cvector_t *u = cvector_alloc(N);
  cvector_t *l = cvector_alloc(N);
  if (!u || !l) {
    cvector_free(u);
    cvector_free(l);

    return;
  }

  // Gaussian * (1, i)/\sqrt(2) at k0 = 0
  dirac_packet_1d(u, l, DX, 0.0, 0.0, 5.0, 0, 1.0, 1.0, 1.0);
  for (int i = 0; i < N; i++) {
    complex_t g = u->data[i];

    u->data[i] = c_scale(g, 1.0 / sqrt(2.0));
    l->data[i] = c_scale(c_new(-g.im, g.re), 1.0 / sqrt(2.0));
  }

  printf("  %8s %10s %10s\n", "t", "<x>", "analytic");
  const double dt = M_PI / 400.0;
  for (int q = 0; q <= 8; q++) {
    double t = q * M_PI / 8.0;

    printf("  %8.4f %10.4f %10.4f\n", t, dirac_position_1d(u, l, DX),
           -0.5 * (1.0 - cos(2.0 * t)));
    dirac_evolve_1d(u, l, NULL, DX, dt, 50, 1.0, 1.0, 1.0);
  }

  cvector_free(u);
  cvector_free(l);
}

static void klein_demo(void) {
  printf("  === Klein step: E = 2 m c^2 ===\n");
  printf("  %6s  %10s  %10s   regime\n", "V0", "T (sim)", "T (exact)");

  const double heights[] = {0.0, 0.5, 2.0, 4.0, 6.0, 10.0, 20.0};
  const int nh = (int)(sizeof heights / sizeof heights[0]);

  for (int h = 0; h < nh; h++) {
    cvector_t *u = cvector_alloc(N);
    cvector_t *l = cvector_alloc(N);
    double *v = malloc(N * sizeof *v);

    if (!u || !l || !v) {
      cvector_free(u);
      cvector_free(l);
      free(v);

      return;
    }

    for (int i = 0; i < N; i++) {
      v[i] = (i > N / 2) ? heights[h] : 0.0;
    }

    dirac_packet_1d(u, l, DX, -50.0, sqrt(3.0), 6.0, 1, 1.0, 1.0, 1.0);
    dirac_evolve_1d(u, l, v, DX, 0.04, 2500, 1.0, 1.0, 1.0);

    double right = 0.0;
    for (int i = N / 2 + 1; i < N; i++) {
      right += c_abs2(u->data[i]) + c_abs2(l->data[i]);
    }

    right *= DX;

    double exact = dirac_step_transmission(2.0, heights[h], 1.0, 1.0, 1.0);
    const char *regime = (heights[h] < 1.0)   ? "normal"
                         : (heights[h] < 3.0) ? "no propagating state"
                                              : "Klein";

    printf("  %6.1f  %10.4f  %10.4f   %s\n", heights[h], right, exact, regime);

    cvector_free(u);
    cvector_free(l);
    free(v);
  }
}

int main(void) {
  printf(" > Dirac wavepacket dynamics\n\n");

  zitterbewegung_demo();
  klein_demo();

  return 0;
}
