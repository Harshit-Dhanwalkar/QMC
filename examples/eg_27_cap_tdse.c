/*
 * Complex Absorbing Potential (CAP) and Time-Dependent V(x,t)
 *
 * A free Gaussian wavepacket travels toward right edge of grid.
 */

#include "../core/complex.h"
#include "../core/ode/crank_nicolson.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double V_free(double x, double t, void *params) { return 0.0; }

int main(void) {
  printf(" > Complex Absorbing Potential (CAP) and Time-Dependent TDSE\n\n");

  int N = 400;
  double x_min = -20.0, x_max = 20.0;
  double dx = (x_max - x_min) / (N - 1);
  double hbar_sq_2m = 0.5;
  double dt = 0.02;
  int total_steps = 400;
  int report_every = 40;

  double *x = malloc(N * sizeof *x);
  for (int i = 0; i < N; i++) {
    x[i] = x_min + i * dx;
  }

  double *absorb = malloc(N * sizeof *absorb);
  cap_build_monomial(x, N, 4.0, 5.0, 2, absorb);

  cvector_t *psi = cvector_alloc(N);
  double sigma0 = 1.0, k0 = 3.0, x0 = -10.0;
  for (int i = 0; i < N; i++) {
    double dxg = x[i] - x0;
    double env = exp(-dxg * dxg / (4.0 * sigma0 * sigma0));

    psi->data[i].re = env * cos(k0 * x[i]);
    psi->data[i].im = env * sin(k0 * x[i]);
  }

  printf("   step    t        norm (expect decay as packet enters CAP)\n");
  printf("   ----   ------   --------\n");

  for (int block = 0; block < total_steps / report_every; block++) {
    crank_nicolson_evolve_time_dependent(x, N, dx, hbar_sq_2m, V_free, NULL,
                                         absorb, psi, block * report_every * dt,
                                         dt, report_every);

    double n = 0.0;
    for (int i = 0; i < N; i++) {
      n += c_abs2(psi->data[i]);
    }
    n *= dx;

    printf("   %4d   %.3f   %.6f\n", (block + 1) * report_every,
           (block + 1) * report_every * dt, n);
  }

  free(x);
  free(absorb);
  cvector_free(psi);

  return 0;
}
