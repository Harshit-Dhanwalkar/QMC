#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/ode/numerov.h"
#include "../core/utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Count nodes
static int count_nodes_robust(const cvector_t *psi, int N, double dx) {
  (void)dx;
  double peak = 0.0;
  for (int i = 0; i < N; i++)
    if (fabs(psi->data[i].re) > peak)
      peak = fabs(psi->data[i].re);
  if (peak < 1e-300)
    return 0;
  double fl = peak * 1e-6;
  // Collect significant values
  double *sig = malloc(N * sizeof *sig);
  int nsig = 0;
  for (int i = 0; i < N; i++)
    if (fabs(psi->data[i].re) > fl)
      sig[nsig++] = psi->data[i].re;
  int nodes = 0;
  for (int i = 1; i < nsig; i++)
    if (sig[i - 1] * sig[i] < 0.0)
      nodes++;
  free(sig);
  return nodes;
}

int main() {
  printf(" > Testing Numerov shooting for harmonic oscillator...\n");

  int N = 201;
  double x_min = -4.0, x_max = 4.0;
  double dx = (x_max - x_min) / (N - 1);

  double *x = linspace(x_min, x_max, N);
  if (!x) {
    printf("FAIL: memory\n");
    return 1;
  }
  double *V = malloc(N * sizeof(double));
  for (int i = 0; i < N; i++)
    V[i] = 0.5 * x[i] * x[i];

  numerov_params_t params = {
      .x = x, .V = V, .n = N, .dx = dx, .hbar_sq_2m = 0.5};

  // Ground state
  numerov_solution_t *sol = numerov_shoot(&params, 0.0, 1e-10);
  if (!sol) {
    printf("FAIL: no ground state solution\n");
    free(x);
    free(V);
    return 1;
  }

  printf("   Found energy: %f (expected 0.5)\n", sol->energy);
  if (fabs(sol->energy - 0.5) > 0.01) {
    printf("FAIL: energy mismatch.\n");
    numerov_solution_free(sol);
    free(x);
    free(V);
    return 1;
  }

  double norm_sq = 0.0, x2 = 0.0;
  for (int i = 0; i < N; i++) {
    norm_sq += sol->psi->data[i].re * sol->psi->data[i].re;
    x2 += x[i] * x[i] * sol->psi->data[i].re * sol->psi->data[i].re;
  }
  norm_sq *= dx;
  x2 *= dx;
  printf("   Normalization: %f\n", norm_sq);

  // Wavefunction vs analytic (in allowed region |x|<1.5)
  double pi_q = pow(M_PI, -0.25), max_err = 0.0;
  int nchk = 0;
  for (int i = 0; i < N; i++) {
    if (fabs(x[i]) > 1.5)
      continue;
    double ana = pi_q * exp(-x[i] * x[i] / 2.0);
    if (fabs(ana) > 1e-3) {
      double err = fabs(sol->psi->data[i].re - ana) / fabs(ana);
      if (err > max_err)
        max_err = err;
      nchk++;
    }
  }
  printf("   Max rel error vs analytic (|x|<1.5): %.2e (%d pts)\n", max_err,
         nchk);
  printf("   <x^2>: %f (expected 0.5)\n", x2);
  numerov_solution_free(sol);

  // First excited state
  sol = numerov_shoot(&params, 1.0, 1e-10);
  if (!sol) {
    printf("FAIL: no n=1 solution\n");
    free(x);
    free(V);
    return 1;
  }
  printf("   E_1 = %f (expected 1.5)\n", sol->energy);
  int nodes = count_nodes_robust(sol->psi, N, dx);
  printf("   n=1 nodes: %d (expected 1)\n", nodes);
  int pass = (fabs(sol->energy - 0.5) < 0.01 || 1) // energy already checked
             && fabs(sol->energy - 1.5) < 0.01 && nodes == 1 &&
             fabs(norm_sq - 1.0) < 0.01 && fabs(x2 - 0.5) < 0.02 &&
             max_err < 0.02;
  numerov_solution_free(sol);
  free(x);
  free(V);

  if (pass) {
    printf("   Numerov test passed.\n");
    return 0;
  }
  printf("FAIL\n");
  return 1;
}
