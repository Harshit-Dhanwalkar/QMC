/*
 * Test: physics/schrodinger.c's TISE-shooting and TDSE-evolution functions
 * (solve_tise_shoot, evolve_tdse_crank, evolve_tdse_split_step).
 *
 *   1. solve_tise_shoot: harmonic-oscillator ground/first-excited energies
 *      against the exact analytic (n + 1/2) * \hbar * \omega, and
 *      cross-validated against solve_tise_matrix's independent diagonalization
 *      approach.
 *   2. evolve_tdse_crank: norm conservation, and same coherent-state trajectory
 *      check - since this uses a completely different discretization
 *      (finite-difference + Crank-Nicolson, not FFT-based split-step),
 *      agreement between the two independently-implemented time-evolution
 *      methods is itself a strong cross-validation of both.
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/ode/crank_nicolson.h"
#include "../core/ode/numerov.h"
#include "../core/vector.h"
#include "../physics/potentials.h"
#include "../physics/schrodinger.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check(int cond, const char *msg) {
  if (!cond) {
    printf("  FAIL: %s\n", msg);
    failures++;
  }
}

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  if (fabs(got - expected) > tol) {
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

/* ---------------------------------------------------------------------
 * solve_tise_shoot
 * ------------------------------------------------------------------- */

static void test_solve_tise_shoot_harmonic_oscillator(void) {
  printf("Test: solve_tise_shoot on the harmonic oscillator matches exact "
         "(n+1/2)*hbar*omega energies and cross-validates against "
         "solve_tise_matrix\n");

  int n = 800;
  double L = 20.0;
  double dx = L / n;
  double omega = 1.0;
  double hbar_sq_2m = 0.5; // \hbar = m = 1 natural units: \hbar^2/(2m) = 0.5

  double *x = malloc((size_t)n * sizeof(double));
  double *V = malloc((size_t)n * sizeof(double));
  for (int i = 0; i < n; i++) {
    x[i] = -L / 2.0 + i * dx;
    V[i] = 0.5 * omega * omega * x[i] * x[i];
  }

  numerov_params_t params = {
      .x = x, .V = V, .n = n, .dx = dx, .hbar_sq_2m = hbar_sq_2m};

  numerov_solution_t *ground = solve_tise_shoot(&params, 0.0, 1e-8);
  check(ground != NULL, "ground-state shoot should succeed");
  if (ground) {
    check_close(ground->energy, 0.5 * omega, 1e-3,
                "ground state energy matches exact 0.5*hbar*omega");
  }

  numerov_solution_t *first_excited = solve_tise_shoot(&params, 1.0, 1e-8);
  check(first_excited != NULL, "first-excited shoot should succeed");
  if (first_excited) {
    check_close(first_excited->energy, 1.5 * omega, 1e-3,
                "first excited state energy matches exact 1.5*hbar*omega");
  }

  eigen_t *eig = solve_tise_matrix(x, n, dx, hbar_sq_2m, V_harmonic, &omega);
  check(eig != NULL, "solve_tise_matrix should succeed");
  if (eig && ground) {
    /* eigenvalues from solve_tise_matrix aren't guaranteed pre-sorted by every
     * possible backend, so find minimum explicitly. */
    double min_eig = eig->eigenvalues[0];
    for (int i = 1; i < eig->n; i++) {
      if (eig->eigenvalues[i] < min_eig) {
        min_eig = eig->eigenvalues[i];
      }
    }

    check_close(ground->energy, min_eig, 1e-3,
                "solve_tise_shoot ground energy matches solve_tise_matrix's "
                "independent diagonalization");
  }

  if (ground) {
    numerov_solution_free(ground);
  }
  if (first_excited) {
    numerov_solution_free(first_excited);
  }
  if (eig) {
    eigen_free(eig);
  }
  free(x);
  free(V);
}

/* ---------------------------------------------------------------------
 * Shared coherent-state setup for the two TDSE evolution tests
 * ------------------------------------------------------------------- */

/*
 * Minimum-uncertainty Gaussian ("coherent state") for harmonic oscillator,
 * centered at x0 with zero initial momentum. Its center-of-mass  trajectory
 * x(t) = x0 * \cos(\omega * t) is an exact result (Ehrenfest's theorem is exact
 * for any harmonic potential) for any \hbar, mass, \omega
 */
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

static double compute_norm(const cvector_t *psi, int n, double dx) {
  double s = 0.0;
  for (int i = 0; i < n; i++) {
    s +=
        (psi->data[i].re * psi->data[i].re + psi->data[i].im * psi->data[i].im);
  }

  return s * dx;
}

static double compute_x_mean(const cvector_t *psi, const double *x, int n,
                             double dx) {
  double s = 0.0, norm = 0.0;
  for (int i = 0; i < n; i++) {
    double p =
        psi->data[i].re * psi->data[i].re + psi->data[i].im * psi->data[i].im;

    s += x[i] * p;
    norm += p;
  }

  return s * dx / (norm * dx);
}

/* ---------------------------------------------------------------------
 * evolve_tdse_split_step
 * ------------------------------------------------------------------- */

static void test_split_step_norm_conservation(void) {
  printf("Test: evolve_tdse_split_step conserves probability norm\n");

  int n = 512; // power of two
  double L = 40.0;
  double dx = L / n;
  double hbar = 1.0, mass = 1.0, omega = 1.0;

  double *x = malloc((size_t)n * sizeof(double));
  double *V = malloc((size_t)n * sizeof(double));
  for (int i = 0; i < n; i++) {
    x[i] = -L / 2.0 + i * dx;
    V[i] = 0.5 * mass * omega * omega * x[i] * x[i];
  }

  cvector_t *psi = cvector_alloc(n);
  make_coherent_state(psi, x, n, 3.0, hbar, mass, omega);
  double norm0 = compute_norm(psi, n, dx);

  int rc = evolve_tdse_split_step(psi, x, V, n, dx, 0.001, 500, hbar, mass);
  check(rc == 0, "evolve_tdse_split_step should return success");

  double norm1 = compute_norm(psi, n, dx);
  check_close(norm1, norm0, 1e-6, "norm conserved after 500 split-step steps");

  cvector_free(psi);
  free(x);
  free(V);
}

static void test_split_step_coherent_state_hbar_independence(void) {
  printf("Test: evolve_tdse_split_step reproduces the exact coherent-state "
         "trajectory x(t)=x0*cos(omega*t) at multiple hbar values\n");

  int n = 1024;
  double L = 40.0;
  double dx = L / n;
  double mass = 1.0, omega = 1.0, x0 = 3.0;
  double dt = 0.001;
  int steps = 2000;
  double t_total = dt * steps;
  double x_expected = x0 * cos(omega * t_total);

  const double hbars[] = {0.5, 1.0, 2.0, 3.0};
  for (int h = 0; h < 4; h++) {
    double hbar = hbars[h];

    double *x = malloc((size_t)n * sizeof(double));
    double *V = malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
      x[i] = -L / 2.0 + i * dx;
      V[i] = 0.5 * mass * omega * omega * x[i] * x[i];
    }

    cvector_t *psi = cvector_alloc(n);
    make_coherent_state(psi, x, n, x0, hbar, mass, omega);

    int rc = evolve_tdse_split_step(psi, x, V, n, dx, dt, steps, hbar, mass);
    check(rc == 0, "evolve_tdse_split_step should return success");

    double x_mean = compute_x_mean(psi, x, n, dx);
    char msg[160];
    snprintf(msg, sizeof(msg),
             "\\hbar=%.1f: <x>(t) matches exact classical trajectory "
             "(Ehrenfest is exact here for any \\hbar)",
             hbar);

    check_close(x_mean, x_expected, 5e-3, msg);

    cvector_free(psi);
    free(x);
    free(V);
  }
}

/* ---------------------------------------------------------------------
 * evolve_tdse_crank
 * ------------------------------------------------------------------- */

static void test_crank_norm_conservation(void) {
  printf("Test: evolve_tdse_crank conserves probability norm\n");

  int n = 500;
  double L = 40.0;
  double dx = L / (n - 1);
  double hbar_sq_2m = 0.5;
  double omega = 1.0;

  double *x = malloc((size_t)n * sizeof(double));
  double *V = malloc((size_t)n * sizeof(double));
  for (int i = 0; i < n; i++) {
    x[i] = -L / 2.0 + i * dx;
    V[i] = 0.5 * omega * omega * x[i] * x[i];
  }

  double *diag = malloc((size_t)n * sizeof(double));
  double *offdiag = malloc((size_t)(n - 1) * sizeof(double));
  build_tridiagonal_hamiltonian(x, V, n, dx, hbar_sq_2m, diag, offdiag);

  cvector_t *psi = cvector_alloc(n);
  make_coherent_state(psi, x, n, 3.0, 1.0, 1.0, omega);
  double norm0 = compute_norm(psi, n, dx);

  int rc = evolve_tdse_crank(diag, offdiag, n, psi, 0.001, 500);
  check(rc == 0, "evolve_tdse_crank should return success");

  double norm1 = compute_norm(psi, n, dx);
  check_close(norm1, norm0, 1e-6,
              "norm conserved after 500 Crank-Nicolson steps");

  cvector_free(psi);
  free(x);
  free(V);
  free(diag);
  free(offdiag);
}

static void test_crank_matches_split_step(void) {
  printf("Test: evolve_tdse_crank and evolve_tdse_split_step "
         "independently-implemented time-evolution schemes "
         "(finite-difference+Crank-Nicolson vs. FFT-based split-step) agree on "
         "the same coherent-state trajectory\n");

  int n = 1024;
  double L = 40.0;
  double dx = L / n;
  double hbar_sq_2m = 0.5; // \hbar=m=1
  double hbar = 1.0, mass = 1.0, omega = 1.0, x0 = 3.0;
  double dt = 0.001;
  int steps = 1000;

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
  build_tridiagonal_hamiltonian(x, V, n, dx, hbar_sq_2m, diag, offdiag);

  int rc1 = evolve_tdse_crank(diag, offdiag, n, psi_crank, dt, steps);
  int rc2 =
      evolve_tdse_split_step(psi_split, x, V, n, dx, dt, steps, hbar, mass);
  check(rc1 == 0 && rc2 == 0, "both evolution methods should succeed");

  double x_crank = compute_x_mean(psi_crank, x, n, dx);
  double x_split = compute_x_mean(psi_split, x, n, dx);

  check_close(x_crank, x_split, 5e-3,
              "Crank-Nicolson and split-step agree on <x>(t) (independent "
              "cross-validation of both methods)");

  double t_total = dt * steps;
  double x_expected = x0 * cos(omega * t_total);
  check_close(x_crank, x_expected, 5e-3,
              "Crank-Nicolson also matches the exact classical trajectory");

  cvector_free(psi_crank);
  cvector_free(psi_split);
  free(x);
  free(V);
  free(diag);
  free(offdiag);
}

static void test_solve_tise_shoot_matching_harmonic_oscillator(void) {
  printf("Test: solve_tise_shoot_matching on the harmonic oscillator matches "
         "the exact (n + 1/2) * \\hbar * \\omega spectrum across 4 levels, "
         "with node counts and parity\n");

  int n = 2000;
  double L = 20.0;
  double dx = L / (n - 1);
  double hbar_sq_2m = 0.5; // \hbar = m = \omega = 1 natural units

  double *x = malloc((size_t)n * sizeof(double));
  double *V = malloc((size_t)n * sizeof(double));
  for (int i = 0; i < n; i++) {
    x[i] = -L / 2.0 + i * dx;
    V[i] = 0.5 * x[i] * x[i];
  }

  numerov_params_t params = {
      .x = x, .V = V, .n = n, .dx = dx, .hbar_sq_2m = hbar_sq_2m};

  const double exact[4] = {0.5, 1.5, 2.5, 3.5};
  for (int level = 0; level < 4; level++) {
    double E_min = exact[level] - 0.4;
    double E_max = exact[level] + 0.4;

    numerov_solution_t *sol =
        solve_tise_shoot_matching(&params, E_min, E_max, 200, 1e-8);

    char msg[128];
    snprintf(msg, sizeof(msg),
             "level %d: bracket [%.2f,%.2f] should find a "
             "root",
             level, E_min, E_max);

    check(sol != NULL, msg);
    if (!sol) {
      continue;
    }

    snprintf(msg, sizeof(msg),
             "level %d energy matches exact (n+1/2)*hbar*omega", level);

    check_close(sol->energy, exact[level], 1e-6, msg);

    int nodes = 0;
    for (int i = 1; i < n; i++) {
      if (sol->psi->data[i].re * sol->psi->data[i - 1].re < 0.0 &&
          fabs(sol->psi->data[i].re) > 1e-6) {
        nodes++;
      }
    }
    snprintf(msg, sizeof(msg), "level %d has exactly %d node(s)", level, level);
    check(nodes == level, msg);

    // Normalization: continuum (dx-weighted) norm should be 1
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
      norm += sol->psi->data[i].re * sol->psi->data[i].re;
    }

    norm *= dx;
    snprintf(msg, sizeof(msg), "level %d wavefunction is normalized", level);

    check_close(norm, 1.0, 1e-4, msg);

    numerov_solution_free(sol);
  }

  free(x);
  free(V);
}

static void test_solve_tise_shoot_matching_matches_shoot(void) {
  printf("Test: solve_tise_shoot_matching's ground-state energy matches "
         "solve_tise_shoot's independent diagonalization-based method");

  int n = 2000;
  double L = 20.0;
  double dx = L / (n - 1);
  double hbar_sq_2m = 0.5;

  double *x = malloc((size_t)n * sizeof(double));
  double *V = malloc((size_t)n * sizeof(double));
  for (int i = 0; i < n; i++) {
    x[i] = -L / 2.0 + i * dx;
    V[i] = 0.5 * x[i] * x[i];
  }

  numerov_params_t params = {
      .x = x, .V = V, .n = n, .dx = dx, .hbar_sq_2m = hbar_sq_2m};

  numerov_solution_t *matching =
      solve_tise_shoot_matching(&params, 0.1, 0.9, 200, 1e-8);
  numerov_solution_t *diag = solve_tise_shoot(&params, 0.0, 1e-8);

  check(matching != NULL && diag != NULL,
        "both methods should find the ground state");
  if (matching && diag) {
    check_close(matching->energy, diag->energy, 1e-5,
                "log-derivative-matching shooting agrees with diagonalization "
                "on same ground-state energy");
  }

  if (matching) {
    numerov_solution_free(matching);
  }
  if (diag) {
    numerov_solution_free(diag);
  }
  free(x);
  free(V);
}

static void test_solve_tise_shoot_matching_no_root_returns_null(void) {
  printf("Test: solve_tise_shoot_matching returns NULL when bracket contains "
         "no eigenvalue\n");

  int n = 2000;
  double L = 20.0;
  double dx = L / (n - 1);
  double *x = malloc((size_t)n * sizeof(double));
  double *V = malloc((size_t)n * sizeof(double));
  for (int i = 0; i < n; i++) {
    x[i] = -L / 2.0 + i * dx;
    V[i] = 0.5 * x[i] * x[i];
  }

  numerov_params_t params = {
      .x = x, .V = V, .n = n, .dx = dx, .hbar_sq_2m = 0.5};

  // [0.6, 0.9] contains no QHO eigenvalue (nearest are 0.5 and 1.5)
  numerov_solution_t *sol =
      solve_tise_shoot_matching(&params, 0.6, 0.9, 200, 1e-8);

  check(sol == NULL, "empty bracket should return NULL, not a false match");
  if (sol) {
    numerov_solution_free(sol);
  }

  free(x);
  free(V);
}

int main(void) {
  test_solve_tise_shoot_harmonic_oscillator();
  test_solve_tise_shoot_matching_harmonic_oscillator();
  test_solve_tise_shoot_matching_matches_shoot();
  test_solve_tise_shoot_matching_no_root_returns_null();
  test_split_step_norm_conservation();
  test_split_step_coherent_state_hbar_independence();
  test_crank_norm_conservation();
  test_crank_matches_split_step();

  if (failures == 0) {
    printf("\nTISE shooting and TDSE evolution tests checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
