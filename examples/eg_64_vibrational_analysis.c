/*
 * Harmonic Vibrational Frequency Analysis
 *
 * Gram-Schmidt on the mass-weighted translation/rotation vectors), and reports
 * the vibrational frequencies in wavenumbers.
 *
 * H2/STO-3G: a linear diatomic has only 1 vibrational mode (3*2 degrees of
 * freedom - 3 translation - 2 rotation, third rotation about bond axis having
 * zero moment of inertia).
 */

#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/vibrational.h"
#include <stdio.h>
#include <stdlib.h>

static molecular_system_t *build_h2(const double geom[][3], int n_atoms,
                                    void *userdata) {
  (void)userdata;
  if (n_atoms != 2) {
    return NULL;
  }

  molecular_system_t *sys = malloc(sizeof *sys);
  sys->n_basis = 2;
  sys->basis = malloc(2 * sizeof *sys->basis);
  sys->basis[0] = molint_basis_sto3g_h(geom[0]);
  sys->basis[1] = molint_basis_sto3g_h(geom[1]);

  const double charges[2] = {1.0, 1.0};
  double centers[2][3];
  for (int a = 0; a < 2; a++) {
    for (int d = 0; d < 3; d++) {
      centers[a][d] = geom[a][d];
    }
  }

  sys->mol = molecule_alloc(2, charges, centers);
  sys->atom_of_basis = malloc(2 * sizeof *sys->atom_of_basis);
  sys->atom_of_basis[0] = 0;
  sys->atom_of_basis[1] = 1;
  sys->n_electrons = 2;

  return sys;
}

static void free_h2(molecular_system_t *sys, void *userdata) {
  (void)userdata;

  basis_function_free(sys->basis[0]);
  basis_function_free(sys->basis[1]);
  free(sys->basis);
  free(sys->atom_of_basis);
  molecule_free(sys->mol);
  free(sys);
}

int main(void) {
  printf(" > Harmonic Vibrational Frequency Analysis: H2/STO-3G\n\n");

  double R = 1.4; // bohr, close to the STO-3G/RHF equilibrium bond length
  double geom0[2][3] = {{0, 0, 0}, {0, 0, R}};

  printf("  Step 1: numerical Hessian (central difference of the analytic "
         "RHF gradient)\n\n");
  double *H = molecular_rhf_hessian_numerical(build_h2, free_h2, NULL, geom0, 2,
                                              1e-12, 300, 1e-3);
  if (!H) {
    fprintf(stderr, "  Hessian computation failed.\n");

    return 1;
  }

  printf("  6x6 Hessian (Hartree/bohr^2):\n");
  for (int i = 0; i < 6; i++) {
    printf("   ");
    for (int j = 0; j < 6; j++) {
      printf(" %+8.4f", H[i * 6 + j]);
    }

    printf("\n");
  }

  printf("\n");

  printf("  Step 2: mass-weight, project out translation/rotation, "
         "diagonalize\n\n");
  double masses[2] = {1.00782503207, 1.00782503207}; // precise H-1 mass

  molecular_vib_result_t *vib =
      molecular_vibrational_analysis(H, 2, masses, geom0);
  free(H);

  if (!vib) {
    fprintf(stderr, "  Vibrational analysis failed.\n");

    return 1;
  }

  printf("  NOTE: %d translation/rotation modes projected out (3 translation + "
         "2 rotation for this linear diatomic - third rotation, about bond "
         "axis, has zero moment of inertia and drops out automatically)\n",
         vib->n_trans_rot);

  printf("  %d vibrational mode(s):\n", vib->n_modes);
  for (int i = 0; i < vib->n_modes; i++) {
    printf("    mode %d: %.2f cm^-1\n", i, vib->frequencies_cm1[i]);
  }

  molecular_vib_result_free(vib);

  return 0;
}
