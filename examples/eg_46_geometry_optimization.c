/*
 * Analytic RHF Nuclear Gradient / Forces (Pulay 1969) : Geometry Optimization
 *
 * NOTE: Every molecular_hf.c calculation has needed a bond length supplied by
 * hand (H2 at the theoretical R=1.4 bohr, LiH at its experimental R=3.015
 * bohr). This example uses molecular_rhf_gradient() an analytic derivative of
 * RHF energy with respect to nuclear coordinates, built from the
 * McMurchie-Davidson "increment/decrement" identity (Reference: Helgaker,
 * Jorgensen & Olsen eq. 9.3.8: differentiating a GTO wrt its own center is
 * itself expressible via the same integral routines with angular momentum
 * shifted by +-1) plus the Pulay energy-weighted-density-matrix term for
 * AO-basis-non-orthogonality (overlap) piece to instead find the equilibrium
 * bond length directly, by walking downhill on the force.
 *
 * Two systems: H2 (for comparison point) and LiH (heteronuclear, p-orbitals,
 * multiple basis functions sharing a center)
 */

#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Returns z-component force on atom 1 (bond-length direction) at bond length R,
 * and optionally the total energy. */
static double h2_force(double R, double *E_out) {
  double c0[3] = {0, 0, 0}, c1[3] = {0, 0, R};

  basis_function_t *h0 = molint_basis_sto3g_h(c0);
  basis_function_t *h1 = molint_basis_sto3g_h(c1);
  basis_function_t *basis[2] = {h0, h1};

  const double charges[2] = {1.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 2, mol, 2, 1e-12, 200);
  if (E_out) {
    *E_out = res->total_energy;
  }
  const int atom_of_basis[2] = {0, 1};
  double *grad = molecular_rhf_gradient(basis, 2, mol, atom_of_basis, res);
  double force_z = -grad[5]; // force = -gradient

  free(grad);
  molecular_hf_result_free(res);
  molecule_free(mol);
  basis_function_free(h0);
  basis_function_free(h1);

  return force_z;
}

static double lih_force(double R, double *E_out) {
  double cLi[3] = {0, 0, 0}, cH[3] = {0, 0, R};

  basis_function_t *li[5];
  molint_basis_sto3g_li(cLi, li);
  basis_function_t *h = molint_basis_sto3g_h(cH);
  basis_function_t *basis[6] = {li[0], li[1], li[2], li[3], li[4], h};

  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);
  molecular_hf_result_t *res = molecular_rhf(basis, 6, mol, 4, 1e-12, 300);
  if (E_out) {
    *E_out = res->total_energy;
  }
  const int atom_of_basis[6] = {0, 0, 0, 0, 0, 1};
  double *grad = molecular_rhf_gradient(basis, 6, mol, atom_of_basis, res);
  double force_z = -grad[5];

  free(grad);
  molecular_hf_result_free(res);
  molecule_free(mol);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }
  basis_function_free(h);

  return force_z;
}

/* Simple steepest-descent-on-the-force optimizer: at equilibrium the force
 * is zero, so a small step along the force direction with a fixed step size
 * converges to the minimum for these smooth, convex-near-the-minimum 1D
 * potential energy curves. */
static double optimize_bond_length(double (*force_fn)(double, double *),
                                   double R0, double step, int max_iter,
                                   double tol, const char *name) {
  double R = R0;
  printf("  Optimizing %s bond length (steepest descent on force):\n", name);

  for (int it = 0; it < max_iter; it++) {
    double E;
    double F = force_fn(R, &E);

    printf("    iter %2d: R=%.6f bohr  E=%.8f Hartree  F=%+.6f\n", it, R, E, F);
    if (fabs(F) < tol) {
      printf("    converged (|F| < %.1e)\n", tol);

      break;
    }

    R += step * F; /* move along the force (downhill in energy) */
  }

  return R;
}

int main(void) {
  printf("=== Analytic RHF Nuclear Gradient : Geometry Optimization ===\n\n");

  double R_h2 = optimize_bond_length(h2_force, 1.0, 0.5, 30, 1e-6, "H2/STO-3G");
  printf("\n  H2/STO-3G equilibrium bond length: %.4f bohr (%.4f Angstrom)\n",
         R_h2, R_h2 * 0.529177);
  printf("  (Literature STO-3G/RHF equilibrium is ~1.34-1.40 bohr)\n\n");

  double R_lih =
      optimize_bond_length(lih_force, 2.5, 1.0, 100, 1e-6, "LiH/STO-3G");
  printf("\n  LiH/STO-3G equilibrium bond length: %.4f bohr (%.4f "
         "Angstrom)\n",
         R_lih, R_lih * 0.529177);
  printf("  (Experimental LiH equilibrium is R=3.015 bohr / 1.596 "
         "Angstrom; STO-3G/RHF\n"
         "   is a small, non-polarized basis so some deviation from "
         "experiment is expected)\n");

  return 0;
}
