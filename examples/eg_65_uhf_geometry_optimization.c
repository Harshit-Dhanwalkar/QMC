/*
 * Analytic UHF Nuclear Gradient: Open-Shell Geometry Optimization
 *
 * Target: LiH+ (the LiH cation, 3 electrons, doublet - 2 \alpha, 1 \beta).
 */

#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double lih_cation_force(double R, double *E_out) {
  double cLi[3] = {0, 0, 0}, cH[3] = {0, 0, R};

  basis_function_t *li[5];
  molint_basis_sto3g_li(cLi, li);
  basis_function_t *h = molint_basis_sto3g_h(cH);
  basis_function_t *basis[6] = {li[0], li[1], li[2], li[3], li[4], h};

  const double charges[2] = {3.0, 1.0};
  double centers[2][3] = {{0, 0, 0}, {0, 0, R}};
  molecule_t *mol = molecule_alloc(2, charges, centers);
  // LiH+: 3 electrons, doublet (2 \alpha, 1 \beta)
  molecular_uhf_result_t *res = molecular_uhf(basis, 6, mol, 2, 1, 1e-12, 300);
  if (E_out) {
    *E_out = res->total_energy;
  }
  const int atom_of_basis[6] = {0, 0, 0, 0, 0, 1};
  double *grad = molecular_uhf_gradient(basis, 6, mol, atom_of_basis, res);
  double force_z = -grad[5]; // force = -gradient

  free(grad);
  molecular_uhf_result_free(res);
  molecule_free(mol);
  for (int i = 0; i < 5; i++) {
    basis_function_free(li[i]);
  }
  basis_function_free(h);

  return force_z;
}

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

    R += step * F;
  }

  return R;
}

int main(void) {
  printf(" > Analytic UHF Nuclear Gradient: Open-Shell Geometry Optimization "
         "\n\n");

  double R = optimize_bond_length(lih_cation_force, 2.5, 3.0, 150, 1e-6,
                                  "LiH+/STO-3G (doublet)");
  printf("\n  LiH+/STO-3G equilibrium bond length: %.4f bohr (%.4f Angstrom)\n",
         R, R * 0.529177);
  printf(
      "  Removing an electron from LiH's closed-shell 2-\\sigma bonding "
      "orbital weakens the bond, so LiH+'s equilibrium bond length is expected "
      "to sit longer than LiH's own ~3.0-3.5 bohr STO-3G/RHF equilibrium (see "
      "eg_47_geometry_optimization.c) - less bonding charge density between "
      "nuclei means less pull holding them together.\n");

  return 0;
}
