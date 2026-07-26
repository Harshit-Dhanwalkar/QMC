/*
Test: klein_gordon_1d_self_consistent (physics/relativistic.c)

1. Constant V(x): self-consistent single-level solver must agree with fast
   klein_gordon_1d (V_avg) solver's corresponding level, since for constant V
   both are solving identical linear problem.
2. Spatially-varying V(x): verify converged (E,psi) satisfies discretized
   nonlinear Klein-Gordon operator equation,
     -\hbar^2 * c^2 (\psi'' / dx^2) + [m^2 * c^4 + 2 * E * V(x) - V(x)^2] \psi
       = E^2 \psi
   directly, via residual norm (independent of how the solver got there).
*/

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/relativistic.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

static int test_constant_V_matches_fast_solver(void) {
  int N = 200;
  double x_min = -5.0, x_max = 5.0;
  double dx = (x_max - x_min) / (N - 1);
  double m = 1.0, hbar = 1.0, c = 1.0;
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
    fail |= check_close(sc->energy, E_ground_fast, 1e-8,
                        "E (self-consistent vs fast, constant V)");
    fail |= !sc->converged;
    klein_gordon_solution_free(sc);
  }

  eigen_free(fast);
  free(x);
  free(V);

  return fail;
}

static int test_spatially_varying_residual(void) {
  int N = 300;
  double x_min = -10.0, x_max = 10.0;
  double dx = (x_max - x_min) / (N - 1);
  double m = 1.0, hbar = 1.0, c = 1.0;

  double *x = malloc(N * sizeof *x);
  double *V = malloc(N * sizeof *V);
  for (int i = 0; i < N; i++) {
    x[i] = x_min + i * dx;
    // Gaussian potential bump
    V[i] = 0.5 * exp(-x[i] * x[i] / 8.0);
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

int main(void) {
  int failed = 0;

  printf("Constant V(x): self-consistent solver matches fast solver:\n");
  failed += test_constant_V_matches_fast_solver();

  printf("Spatially-varying V(x): operator-equation residual check:\n");
  failed += test_spatially_varying_residual();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
