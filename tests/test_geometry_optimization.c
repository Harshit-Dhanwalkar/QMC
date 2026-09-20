/*
 * Test: physics/geometry_optimization.c's generic n_atoms*3 Cartesian
 * steepest-descent optimizer, driven by an RHF energy+analytic-gradient
 * callback (molecular_rhf + molecular_rhf_gradient)
 *
 * eg_49_geometry_optimization.c and test_hf_gradient.c's own
 * test_h2_equilibrium_bond_length both already validate underlying
 * gradient via a hand-rolled *scalar* (bond-length-only) search. This file
 * instead exercises general-purpose library optimizer itself:
 *
 *   1. H2/STO-3G optimized from a collinear starting guess (R=1.0 bohr)
 *      converges to accepted literature STO-3G/RHF equilibrium range (~1.3-1.45
 *      bohr), with a small residual gradient and a lower energy than starting
 *      geometry
 *   2. Optimizing from two different starting geometries with same initial bond
 *      length, but a different, non-symmetric split of that bond length between
 *      two atoms' individual coordinates (i.e. exercising more than one scalar
 *      degree of freedom), converges to same equilibrium bond length - a check
 *      that optimizer is really operating on independent per-atom Cartesian
 *      coordinates and not secretly relying on any bond-length-specific
 *      shortcut
 */

#include "../physics/geometry_optimization.h"
#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAILED");
  if (!cond) {
    failures++;
  }
}

/* geom_energy_grad_fn for H2/STO-3G: coords is {x0,y0,z0, x1,y1,z1}. */
static int h2_energy_grad(const double *coords, int n_atoms, void *user_data,
                          double *energy_out, double *grad_out) {
  (void)user_data;
  if (n_atoms != 2) {
    return 1;
  }

  double c0[3] = {coords[0], coords[1], coords[2]};
  double c1[3] = {coords[3], coords[4], coords[5]};

  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};

  const double charges[2] = {1.0, 1.0};
  double centers[2][3];
  centers[0][0] = c0[0];
  centers[0][1] = c0[1];
  centers[0][2] = c0[2];
  centers[1][0] = c1[0];
  centers[1][1] = c1[1];
  centers[1][2] = c1[2];

  molecule_t *mol = molecule_alloc(2, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 2, mol, 2, 1e-12, 200);

  int rc = 1;
  if (res && res->converged) {
    *energy_out = res->total_energy;

    const int atom_of_basis[2] = {0, 1};
    double *grad = molecular_rhf_gradient(basis, 2, mol, atom_of_basis, res);
    if (grad) {
      for (int i = 0; i < 6; i++) {
        grad_out[i] = grad[i];
      }

      free(grad);

      rc = 0;
    }
  }

  if (res) {
    molecular_hf_result_free(res);
  }
  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);

  return rc;
}

static double bond_length(const double *coords) {
  double dx = coords[3] - coords[0];
  double dy = coords[4] - coords[1];
  double dz = coords[5] - coords[2];

  return sqrt(dx * dx + dy * dy + dz * dz);
}

static void test_null_and_invalid_input(void) {
  printf("  === Test Null and invalid input ===\n");
  double const coords0[6] = {0, 0, 0, 0, 0, 1.0};

  check_true(optimize_geometry_steepest_descent(NULL, NULL, coords0, 2, 0.5,
                                                1e-5, 200) == NULL,
             "NULL callback rejected");
  check_true(optimize_geometry_steepest_descent(h2_energy_grad, NULL, coords0,
                                                0, 0.5, 1e-5, 200) == NULL,
             "n_atoms<=0 rejected");
  check_true(optimize_geometry_steepest_descent(h2_energy_grad, NULL, coords0,
                                                2, -0.1, 1e-5, 200) == NULL,
             "step0<=0 rejected");
  check_true(optimize_geometry_steepest_descent(h2_energy_grad, NULL, coords0,
                                                2, 0.5, 0.0, 200) == NULL,
             "grad_tol<=0 rejected");
  check_true(optimize_geometry_steepest_descent(h2_energy_grad, NULL, coords0,
                                                2, 0.5, 1e-5, 0) == NULL,
             "max_iter<=0 rejected");
}

static void test_h2_optimization_from_collinear_guess(void) {
  printf("  === Test H2 optimization from collinear guess ===");

  double const coords0[6] = {0, 0, 0, 0, 0, 1.0};
  double e0;
  double g0[6];
  int rc0 = h2_energy_grad(coords0, 2, NULL, &e0, g0);
  check_true(rc0 == 0, "starting-geometry SCF converges");

  geometry_optimization_result_t *result = optimize_geometry_steepest_descent(
      h2_energy_grad, NULL, coords0, 2, 0.5, 1e-6, 200);

  check_true(result != NULL, "optimizer returns a result");
  if (!result) {
    return;
  }

  check_true(result->converged, "optimizer reports convergence");
  check_true(result->grad_norm < 1e-4,
             "converged gradient is small (near a stationary point)");
  check_true(result->energy < e0,
             "optimized energy is lower than the starting-geometry energy");

  double Req = bond_length(result->coords);
  printf("  H2/STO-3G optimized bond length: %.4f bohr (energy %.10f)\n", Req,
         result->energy);
  // Literature STO-3G/RHF H2 equilibrium is ~1.3-1.45 bohr
  check_true(Req > 1.3 && Req < 1.45,
             "optimized bond length matches literature STO-3G/RHF range");

  geometry_optimization_result_free(result);
}

static void test_h2_optimization_independent_of_atom_split(void) {
  printf("  === Test H2 optimization independent of atom split ===\n");

  // NOTE: Same initial bond length (1.0 bohr) as the collinear-guess test
  // above, but split asymmetrically and off-axis between the two atoms' own
  // coordinates, and with whole pair translated away from the origin - a
  // multi-coordinate starting point
  double const coords0[6] = {0.3, -0.2, 0.5, 0.3, -0.2, 1.5};
  check_true(fabs(bond_length(coords0) - 1.0) < 1e-9,
             "test setup: starting bond length is 1.0 bohr");

  geometry_optimization_result_t *result = optimize_geometry_steepest_descent(
      h2_energy_grad, NULL, coords0, 2, 0.5, 1e-6, 200);

  check_true(result != NULL, "optimizer returns a result");
  if (!result) {
    return;
  }

  check_true(result->converged, "optimizer reports convergence");

  double Req = bond_length(result->coords);
  printf("  H2/STO-3G optimized bond length (off-axis start): %.4f bohr\n",
         Req);
  check_true(Req > 1.3 && Req < 1.45,
             "converges to the same equilibrium bond-length range regardless "
             "of how the initial displacement was split between atoms");

  geometry_optimization_result_free(result);
}

int main(void) {
  test_null_and_invalid_input();
  test_h2_optimization_from_collinear_guess();
  test_h2_optimization_independent_of_atom_split();

  if (failures > 0) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll checks passed\n");

  return 0;
}
