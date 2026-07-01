#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/ode/numerov.h"
#include "../core/utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
  printf(" > Testing Numerov shooting for harmonic oscillator...\n");

  int N = 1001; // grid size
  double x_min = -8.0, x_max = 8.0;
  double dx = (x_max - x_min) / (N - 1);
  double *x = linspace(x_min, x_max, N);
  if (!x) {
    printf("FAIL: memory\n");
    return 1;
  }

  // Potential: harmonic oscillator (m=1, omega=1, in atomic units)
  double omega = 1.0;
  double mass = 1.0;
  double *V = malloc(N * sizeof(double));
  for (int i = 0; i < N; i++) {
    V[i] = 0.5 * mass * omega * omega * x[i] * x[i];
  }

  numerov_params_t params;
  params.x = x;
  params.V = V;
  params.n = N;
  params.dx = dx;
  params.hbar_sq_2m = HBAR_2M;

  // For atomic units (\hbar=1, m=1), \hbar_sq_2m = 0.5
  // HACK: just set it to 0.5.
  params.hbar_sq_2m = 0.5; // \hbar^2/(2m) with \hbar=1, m=1 -> 0.5

  // Ground state energy should be 0.5
  // TEST:
  // numerov_solution_t *sol = numerov_shoot(&params, 0.5, 1e-8);
  // numerov_solution_t *sol = numerov_shoot(&params, 1.0, 1e-8);
  numerov_solution_t *sol = numerov_shoot(&params, 0.0, 1e-10);
  if (!sol) {
    printf("FAIL: no solution found.\n");
    free(x);
    free(V);
    return 1;
  }
  // DEBUG:
  if (!sol->psi) {
    printf("FAIL: solution has no wavefunction.\n");
    numerov_solution_free(sol);
    free(x);
    free(V);
    return 1;
  }

  printf("    Found energy: %f (expected 0.5)\n", sol->energy);
  if (fabs(sol->energy - 0.5) > 1e-4) {
    printf("FAIL: energy mismatch.\n");
    numerov_solution_free(sol);
    free(x);
    free(V);
    return 1;
  }

  // Check normalization (should be ~1)
  double norm_sq = 0.0;
  for (int i = 0; i < N; i++)
    norm_sq += c_abs2(sol->psi->data[i]);
  norm_sq *= dx;
  printf("    Normalization: %f\n", norm_sq);

  numerov_solution_free(sol);
  free(x);
  free(V);
  return 0;
}
