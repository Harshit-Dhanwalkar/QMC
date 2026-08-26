/*
 * Solving and Evolving the Schrodinger Equation: 4 Methods, One Physical System
 *
 * physics/schrodinger.c exposes five solvers/evolvers; only solve_tise_matrix
 * had any example before this one. This example demonstrates the other four on
 * the same physical system (the 1D quantum harmonic oscillator), cross-checking
 * each against either an exact analytic result or against each other:
 *
 * 1. solve_tise_shoot       - Numerov shooting for a specific bound-state
 *                                energy level, checked against the exact
 *                                (n + 1/2)*hbar*omega spectrum.
 * 2. solve_tise_shoot_matching - log-derivative-matching shooting (distinct
 *                                from solve_tise_shoot, whose underlying
 *                                numerov_shoot() actually diagonalizes a
 *                                tridiagonal matrix despite its name),
 *                                cross-checked against both the exact spectrum
 *                                and solve_tise_shoot's independent
 *                                diagonalization result.
 * 3. evolve_tdse_crank      - Crank-Nicolson time evolution of a
 *                                state, cross-checked against both the exact
 *                                trajectory and evolve_tdse_crank's
 *                                independently-discretized result.
 * 4. evolve_tdse_split_step - FFT-based split-step evolution of same coherent
 *                                state, cross-checked against both the exact
 *                                trajectory and evolve_tdse_crank's
 *                                independently-discretized result.
 */

#include "../core/complex.h"
#include "../core/ode/crank_nicolson.h"
#include "../core/ode/numerov.h"
#include "../core/vector.h"
#include "../physics/schrodinger.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void make_coherent_state(cvector_t *psi, const double *x, int n,
                                double x0, double hbar, double mass,
                                double omega) {
  double sigma2 = hbar / (2.0 * mass * omega);
  double norm_prefactor = pow(2.0 * M_PI * sigma2, -0.25);

  for (int i = 0; i < n; i++) {
    double dx_ = x[i] - x0;
    double val = norm_prefactor * exp(-dx_ * dx_ / (4.0 * sigma2));

    psi->data[i] = c_real(val);
  }
}

static double compute_x_mean(const cvector_t *psi, const double *x, int n) {
  double s = 0.0, norm = 0.0;

  for (int i = 0; i < n; i++) {
    double p =
        psi->data[i].re * psi->data[i].re + psi->data[i].im * psi->data[i].im;

    s += x[i] * p;
    norm += p;
  }

  return s / norm;
}

int main(void) {
  printf(" > 4 Ways to Solve/Evolve the Harmonic Oscillator\n\n");

  double omega = 1.0, mass = 1.0, hbar = 1.0;

  // 1. solve_tise_shoot: bound-state energies via Numerov shooting
  printf("Step 1: solve_tise_shoot - Numerov shooting for bound states\n\n");

  int n_shoot = 800;
  double L_shoot = 20.0;
  double dx_shoot = L_shoot / n_shoot;
  double *x_shoot = malloc((size_t)n_shoot * sizeof(double));
  double *V_shoot = malloc((size_t)n_shoot * sizeof(double));

  for (int i = 0; i < n_shoot; i++) {
    x_shoot[i] = -L_shoot / 2.0 + i * dx_shoot;
    V_shoot[i] = 0.5 * omega * omega * x_shoot[i] * x_shoot[i];
  }

  numerov_params_t params = {.x = x_shoot,
                             .V = V_shoot,
                             .n = n_shoot,
                             .dx = dx_shoot,
                             .hbar_sq_2m = 0.5};

  for (int level = 0; level < 3; level++) {
    numerov_solution_t *sol = solve_tise_shoot(&params, (double)level, 1e-8);
    if (sol) {
      double exact = (level + 0.5) * omega;
      printf("  n=%d: E_shoot = %.6f  (exact = %.6f, diff = %.2e)\n", level,
             sol->energy, exact, fabs(sol->energy - exact));

      numerov_solution_free(sol);
    }
  }

  free(x_shoot);
  free(V_shoot);

  // 1.5. solve_tise_shoot_matching: log-derivative-matching shooting
  printf("\nStep 1.5: solve_tise_shoot_matching - log-derivative-matching "
         "shooting (distinct from solve_tise_shoot's diagonalization-based "
         "numerov_shoot)\n\n");

  int n_match = 2000;
  double L_match = 20.0;
  double dx_match = L_match / (n_match - 1);
  double *x_match = malloc((size_t)n_match * sizeof(double));
  double *V_match = malloc((size_t)n_match * sizeof(double));

  for (int i = 0; i < n_match; i++) {
    x_match[i] = -L_match / 2.0 + i * dx_match;
    V_match[i] = 0.5 * omega * omega * x_match[i] * x_match[i];
  }

  numerov_params_t match_params = {.x = x_match,
                                   .V = V_match,
                                   .n = n_match,
                                   .dx = dx_match,
                                   .hbar_sq_2m = 0.5};

  const double levels[3] = {0.5, 1.5, 2.5};
  for (int level = 0; level < 3; level++) {
    double E_min = levels[level] - 0.4, E_max = levels[level] + 0.4;

    numerov_solution_t *sol =
        solve_tise_shoot_matching(&match_params, E_min, E_max, 200, 1e-8);

    if (sol) {
      printf("  n=%d: E_matching = %.6f  (exact = %.6f, diff = %.2e)\n", level,
             sol->energy, levels[level], fabs(sol->energy - levels[level]));

      numerov_solution_free(sol);
    }
  }

  free(x_match);
  free(V_match);

  // 2 & 3. Time evolution of a coherent-state wavepacket
  printf("\nStep 2: evolve_tdse_crank and evolve_tdse_split_step - "
         "coherent-state wavepacket dynamics\n\n");

  int n = 1024;
  double L = 40.0, x0 = 3.0, dt = 0.001;
  int steps = 3000;
  double dx = L / n;

  double *x = malloc((size_t)n * sizeof(double));
  double *V = malloc((size_t)n * sizeof(double));
  for (int i = 0; i < n; i++) {
    x[i] = -L / 2.0 + i * dx;
    V[i] = 0.5 * mass * omega * omega * x[i] * x[i];
  }

  cvector_t *psi_crank = cvector_alloc(n);
  cvector_t *psi_split = cvector_alloc(n);

  make_coherent_state(psi_crank, x, n, x0, hbar, mass, omega);
  make_coherent_state(psi_split, x, n, x0, hbar, mass, omega);

  double *diag = malloc((size_t)n * sizeof(double));
  double *offdiag = malloc((size_t)(n - 1) * sizeof(double));

  build_tridiagonal_hamiltonian(x, V, n, dx, 0.5, diag, offdiag);

  evolve_tdse_crank(diag, offdiag, n, psi_crank, dt, steps);
  evolve_tdse_split_step(psi_split, x, V, n, dx, dt, steps, hbar, mass);

  double t_total = dt * steps;
  double x_expected = x0 * cos(omega * t_total);
  double x_crank = compute_x_mean(psi_crank, x, n);
  double x_split = compute_x_mean(psi_split, x, n);

  printf("  After t=%.2f (period=%.2f):\n", t_total, 2.0 * M_PI / omega);
  printf("  <x>(t) exact classical:       %.6f\n", x_expected);
  printf("  <x>(t) Crank-Nicolson:        %.6f  (diff %.2e)\n", x_crank,
         fabs(x_crank - x_expected));
  printf("  <x>(t) split-step (FFT):      %.6f  (diff %.2e)\n\n", x_split,
         fabs(x_split - x_expected));

  cvector_free(psi_crank);
  cvector_free(psi_split);
  free(x);
  free(V);
  free(diag);
  free(offdiag);

  return 0;
}
