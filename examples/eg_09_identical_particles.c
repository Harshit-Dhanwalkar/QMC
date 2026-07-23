/*
 * Identical Particles (Slater Determinant / Permanent)
 *
 * Uses particle-in-a-box eigenstates as simple orbital basis, then
 * demonstrates:
 *  - Fermions (Slater determinant): antisymmetric, vanishes exactly
 *    when two particles share a position (Pauli exclusion) or two
 *    orbitals coincide.
 *  - Bosons (permanent): symmetric, no such exclusion.
 */

#include "../core/complex.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../physics/identical.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Identical Particles: Slater Determinant / Permanent\n\n");

  // Particle-in-a-box orbitals (L=1, analytic): \phi_k(x) = \sin(k * \pi * x)
  int M = 50; // grid points
  double L = 1.0;
  double *x = linspace(0.0, L, M);
  if (!x) {
    return 1;
  }

  int n_orbitals = 3;
  cvector_t *orbitals[3];
  for (int k = 0; k < n_orbitals; k++) {
    orbitals[k] = cvector_alloc(M);

    for (int i = 0; i < M; i++) {
      orbitals[k]->data[i] = c_real(sin((k + 1) * M_PI * x[i] / L));
    }
  }

  printf("   Orbitals: particle-in-a-box states n=1,2,3, grid of %d points\n\n",
         M);

  int idx_diff[2] = {10, 30};
  complex_t det_diff = slater_determinant_value(orbitals, 2, idx_diff);
  printf("   Fermions, orbitals {1,2}, distinct positions: det=%.6f\n",
         det_diff.re);

  int idx_same[2] = {20, 20};
  complex_t det_same = slater_determinant_value(orbitals, 2, idx_same);
  printf("   Fermions, orbitals {1,2}, same position: det=%.2e "
         "(should be 0, Pauli exclusion)\n",
         det_same.re);

  complex_t perm_same = bosonic_permanent_value(orbitals, 2, idx_same);
  printf("   Bosons,   orbitals {1,2}, same position: perm=%.6f "
         "(nonzero - no exclusion for bosons)\n\n",
         perm_same.re);

  int idx3[3] = {8, 22, 41};
  complex_t det3 = slater_determinant_value(orbitals, 3, idx3);
  complex_t perm3 = bosonic_permanent_value(orbitals, 3, idx3);
  printf("   3 particles, orbitals {1,2,3}, distinct positions:\n");
  printf("     determinant (fermions) = %.6f\n", det3.re);
  printf("     permanent   (bosons)   = %.6f\n", perm3.re);

  for (int k = 0; k < n_orbitals; k++) {
    cvector_free(orbitals[k]);
  }

  free(x);

  return 0;
}
