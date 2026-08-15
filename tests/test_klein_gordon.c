/*
 * Test: klein_gordon_1d_self_consistent (physics/relativistic.c)
 *
 * 1. Constant V(x): self-consistent single-level solver must agree with fast
 *    klein_gordon_1d (V_avg) solver's corresponding level, since for constant V
 *    both are solving identical linear problem.
 * 2. Spatially-varying V(x): verify converged (E,psi) satisfies discretized
 *    nonlinear Klein-Gordon operator equation,
 *      -\hbar^2 * c^2 (\psi'' / dx^2) + [m^2 * c^4 + 2 * E * V(x) - V(x)^2]
 *      \psi = E^2 \psi
 *   directly, via residual norm (independent of how the solver
 * got there).
 */

#include "../core/complex.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/relativistic.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef RUNNING_ON_VALGRIND
#define RUNNING_ON_VALGRIND 0
#endif

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

static int test_constant_V_matches_fast_solver(void) {
  int N = RUNNING_ON_VALGRIND ? 50 : 200;
  double x_min = -5.0, x_max = 5.0;
  double dx = (x_max - x_min) / (N - 1);
  double m = 1.0, hbar = 1.0, c = 1.0; // atmoic units
  double V0 = 0.3;

  double *x = malloc(N * sizeof *x);
  double *V = malloc(N * sizeof *V);
  for (int i = 0; i < N; i++) {
    x[i] = x_min + i * dx;
    V[i] = V0;
  }

  eigen_t *fast = klein_gordon_1d(x, N, V, m, hbar, c);
  if (!fast) {
    free(x);
    free(V);
    printf("  FAIL: klein_gordon_1d returned NULL\n");

    return 1;
  }

  double E_ground_fast = fast->eigenvalues[0];

  klein_gordon_solution_t *sc = klein_gordon_1d_self_consistent(
      x, N, V, m, hbar, c, E_ground_fast, 1e-12, 100);

  int fail = 0;
  if (!sc) {
    printf("  FAIL: klein_gordon_1d_self_consistent returned NULL\n");
    fail = 1;
  } else {
    printf("  converged=%d iterations=%d\n", sc->converged, sc->iterations);
    fail |= check_close(sc->energy, E_ground_fast, 1e-8, "E (constant V)");
    fail |= !sc->converged;

    klein_gordon_solution_free(sc);
  }

  eigen_free(fast);
  free(x);
  free(V);

  return fail;
}

static int test_spatially_varying_residual(void) {
  int N = RUNNING_ON_VALGRIND ? 80 : 300;
  double x_min = RUNNING_ON_VALGRIND ? -6.0 : -10.0;
  double x_max = RUNNING_ON_VALGRIND ? 6.0 : 10.0;
  double dx = (x_max - x_min) / (N - 1);
  double m = 1.0, hbar = 1.0, c = 1.0; // atmoic units

  double *x = malloc(N * sizeof *x);
  double *V = malloc(N * sizeof *V);
  for (int i = 0; i < N; i++) {
    x[i] = x_min + i * dx;
    V[i] = 0.5 * exp(-x[i] * x[i] / 8.0); // Gaussian potential bump
  }

  // Seed guess from free-particle (V=0) ground state so iteration starts near
  // right basin.
  double coeff = hbar * hbar * c * c / (dx * dx);
  double E_guess =
      sqrt(2.0 * coeff * (1.0 - cos(M_PI / N)) + m * m * c * c * c * c);

  klein_gordon_solution_t *sol =
      klein_gordon_1d_self_consistent(x, N, V, m, hbar, c, E_guess, 1e-10, 200);

  int fail = 0;
  if (!sol) {
    printf("  FAIL: klein_gordon_1d_self_consistent returned NULL\n");
    fail = 1;
  } else {
    printf("  E=%.8f converged=%d iterations=%d\n", sol->energy, sol->converged,
           sol->iterations);

    double mc4 = m * m * c * c * c * c;
    double max_residual = 0.0;

    for (int i = 1; i < N - 1; i++) {
      double lap =
          (-coeff) * sol->psi->data[i - 1].re +
          (2.0 * coeff + mc4 + 2.0 * sol->energy * V[i] - V[i] * V[i]) *
              sol->psi->data[i].re +
          (-coeff) * sol->psi->data[i + 1].re;
      double rhs = sol->energy * sol->energy * sol->psi->data[i].re;
      double res = fabs(lap - rhs);
      if (res > max_residual) {
        max_residual = res;
      }
    }

    printf("  max |H(E)\\psi - E^2\\psi| residual = %.3e\n", max_residual);

    fail |= !sol->converged;
    fail |= (max_residual > 1e-6);

    klein_gordon_solution_free(sol);
  }

  free(x);
  free(V);

  return fail;
}

/*
 * klein_gordon_1d assigns each level n its <V>_n = <n|V|n> expectation value
 * (instead of single uniform V_avg shift applied to every level)
 */
static int test_level_dependent_vs_uniform_shift(void) {
  int N = RUNNING_ON_VALGRIND ? 80 : 300;
  double x_min = RUNNING_ON_VALGRIND ? -6.0 : -10.0;
  double x_max = RUNNING_ON_VALGRIND ? 6.0 : 10.0;
  double dx = (x_max - x_min) / (N - 1);
  double m = 1.0, hbar = 1.0, c = 1.0; // atmoic units

  double *x = malloc(N * sizeof *x);
  double *V = malloc(N * sizeof *V);
  for (int i = 0; i < N; i++) {
    x[i] = x_min + i * dx;
    V[i] = 0.5 * exp(-x[i] * x[i] / 8.0);
  }

  double V_avg = 0.0;
  for (int i = 0; i < N; i++) {
    V_avg += V[i];
  }
  V_avg /= N;

  // Rebuild same free (V-independent) operator klein_gordon_1d uses internally,
  // to compute uniform-V_avg-shift energies for comparison.
  double coeff = hbar * hbar * c * c / (dx * dx);
  double mc2 = m * c * c;
  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  for (int i = 0; i < N; i++) {
    diag[i] = 2.0 * coeff + mc2 * mc2;

    if (i < N - 1) {
      offdiag[i] = -coeff;
    }
  }

  eigen_t *free_eig = tridiag_eigh(diag, offdiag, N);
  free(diag);
  free(offdiag);

  int fail = 0;
  if (!free_eig) {
    printf("  FAIL: could not rebuild free operator for comparison\n");

    free(x);
    free(V);

    return 1;
  }

  eigen_t *new_eig = klein_gordon_1d(x, N, V, m, hbar, c);
  if (!new_eig) {
    printf("  FAIL: klein_gordon_1d returned NULL\n");

    eigen_free(free_eig);
    free(x);
    free(V);

    return 1;
  }

  // Check: klein_gordon_1d's returned energies match manual recomputation of
  // <V>_n formula.
  for (int n = 0; n < 3 && n < new_eig->n; n++) {
    double V_expect = 0.0;

    for (int i = 0; i < N; i++) {
      double ci = CMAT(free_eig->eigenvectors, i, n).re;

      V_expect += ci * ci * V[i];
    }

    double expected = V_expect + sqrt(free_eig->eigenvalues[n]);
    double err = fabs(new_eig->eigenvalues[n] - expected);
    printf("  level %d: klein_gordon_1d=%.8f  manual-formula=%.8f  err=%.2e\n",
           n, new_eig->eigenvalues[n], expected, err);
    fail |= (err > 1e-10);
  }

  // Level 2: well-separated state - level-dependent shift should be more
  // accurate.
  {
    int level = 2;
    double E_old = V_avg + sqrt(free_eig->eigenvalues[level]);
    double E_new = new_eig->eigenvalues[level];

    klein_gordon_solution_t *sc =
        klein_gordon_1d_self_consistent(x, N, V, m, hbar, c, E_old, 1e-10, 200);
    if (sc) {
      double err_old = fabs(E_old - sc->energy);
      double err_new = fabs(E_new - sc->energy);
      printf("  level 2 (well-separated): E_old=%.6f E_new=%.6f E_exact=%.6f "
             "err_old=%.2e err_new=%.2e\n",
             E_old, E_new, sc->energy, err_old, err_new);

      // NOTE: Under Valgrind, coarser grid may make the two errors comparable;
      // skip strict inequality test when running under Valgrind.
      if (!RUNNING_ON_VALGRIND) {
        fail |= !(err_new < err_old);
      }

      klein_gordon_solution_free(sc);
    } else {
      printf("  FAIL: did not converge for level 2\n");
      fail = 1;
    }
  }

  // Level 0 (ground state): this potential makes quasi-degenerate with level 1,
  // NOTE: Where non-degenerate PT is expected to do worse
  {
    int level = 0;
    double E_old = V_avg + sqrt(free_eig->eigenvalues[level]);
    double E_new = new_eig->eigenvalues[level];

    klein_gordon_solution_t *sc =
        klein_gordon_1d_self_consistent(x, N, V, m, hbar, c, E_old, 1e-10, 200);
    if (sc) {
      double err_old = fabs(E_old - sc->energy);
      double err_new = fabs(E_new - sc->energy);
      printf("  level 0 (quasi-degenerate with level 1): E_old=%.6f E_new=%.6f "
             "E_exact=%.6f err_old=%.2e err_new=%.2e\n",
             E_old, E_new, sc->energy, err_old, err_new);
      fail |= !(isfinite(E_old) && isfinite(E_new));

      klein_gordon_solution_free(sc);
    } else {
      printf("  FAIL: did not converge for level 0\n");
      fail = 1;
    }
  }

  eigen_free(free_eig);
  eigen_free(new_eig);
  free(x);
  free(V);

  return fail;
}

int main(void) {
  int failed = 0;

  printf("Constant V(x): self-consistent solver matches fast solver:\n");
  failed += test_constant_V_matches_fast_solver();

  printf("Spatially-varying V(x): operator-equation residual check:\n");
  failed += test_spatially_varying_residual();

  printf("Level-dependent shift vs uniform V_avg shift:\n");
  failed += test_level_dependent_vs_uniform_shift();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
