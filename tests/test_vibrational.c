/*
 * Test: numerical Hessian (central-difference of the analytic RHF gradient) +
 * harmonic vibrational frequency analysis
 *
 * H2/STO-3G at R=1.4 bohr (close to its equilibrium bond length), checked
 * against an independent reference computed RHF Hessian and
 * hessian.thermo.harmonic_analysis:
 *  RHF energy        :  -1.116714325062551 Hartree
 *  harmonic frequency:  5027.0165517 cm^-1 (single vibrational mode, 5
 *                       trans/rot modes at 0)
 */

#include "../physics/molecular_hf.h"
#include "../physics/molecular_integrals.h"
#include "../physics/vibrational.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(const char *label, double got, double expected,
                        double tol) {
  double err = fabs(got - expected);

  printf("  %s: got=%.6f expected=%.6f err=%.2e\n", label, got, expected, err);
  if (err > tol) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

// H2 geometry-builder callback: rebuilds the 2-basis-function H2/STO-3G
// system at whatever geometry the Hessian driver requests.
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
  if (!sys) {
    return;
  }

  basis_function_free(sys->basis[0]);
  basis_function_free(sys->basis[1]);
  free(sys->basis);
  free(sys->atom_of_basis);
  molecule_free(sys->mol);
  free(sys);
}

static void test_h2_hessian_and_frequency(void) {
  printf("test_h2_hessian_and_frequency:\n");

  double R = 1.4;
  double geom0[2][3] = {{0, 0, 0}, {0, 0, R}};

  double *H = molecular_rhf_hessian_numerical(build_h2, free_h2, NULL, geom0, 2,
                                              1e-12, 300, 1e-3);
  check_true(H != NULL, "Hessian computed successfully");
  if (!H) {
    return;
  }

  // Symmetry sanity check: the Hessian must be symmetric (already enforced
  // by construction, but confirm it wasn't accidentally broken).
  int dim = 6;
  double max_asym = 0.0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      double d = fabs(H[i * dim + j] - H[j * dim + i]);

      if (d > max_asym) {
        max_asym = d;
      }
    }
  }
  check_close("Hessian is symmetric", max_asym, 0.0, 1e-12);

  // Precise H-1 atomic mass (NOT the rounded integer - matters at the
  // ~0.5% level for such a light atom).
  double masses[2] = {1.00782503207, 1.00782503207};

  molecular_vib_result_t *vib =
      molecular_vibrational_analysis(H, 2, masses, geom0);
  check_true(vib != NULL, "vibrational analysis succeeded");
  free(H);

  if (!vib) {
    return;
  }

  check_true(
      vib->n_trans_rot == 5,
      "H2 (linear diatomic) has exactly 5 trans/rot modes projected out (3 "
      "translation + 2 rotation; the third rotation, about bond axis, has zero "
      "moment of inertia and drops out of Gram-Schmidt automatically)");
  check_true(vib->n_modes == 1,
             "exactly 1 genuine vibrational mode remains (3*2 - 5 = 1)");

  if (vib->n_modes == 1) {
    check_close("H2/STO-3G harmonic frequency vs independent reference "
                "(numerical vs analytic Hessian, so not exact)",
                vib->frequencies_cm1[0], 5027.0165517, 15.0);
    check_true(vib->frequencies_cm1[0] > 0,
               "frequency is real (H2 at R=1.4 is near a true minimum, not a "
               "saddle point)");
  }

  molecular_vib_result_free(vib);
}

static void test_invalid_input(void) {
  printf("test_invalid_input:\n");

  double geom0[2][3] = {{0, 0, 0}, {0, 0, 1.4}};

  check_true(molecular_rhf_hessian_numerical(NULL, free_h2, NULL, geom0, 2,
                                             1e-10, 100, 1e-3) == NULL,
             "NULL build_system callback rejected");
  check_true(molecular_rhf_hessian_numerical(build_h2, free_h2, NULL, geom0, 0,
                                             1e-10, 100, 1e-3) == NULL,
             "n_atoms<=0 rejected");
  check_true(molecular_rhf_hessian_numerical(build_h2, free_h2, NULL, geom0, 2,
                                             1e-10, 100, -1.0) == NULL,
             "h<=0 rejected");

  check_true(molecular_vibrational_analysis(NULL, 2, NULL, geom0) == NULL,
             "NULL hessian rejected");

  const double dummy[36] = {0};
  check_true(molecular_vibrational_analysis(dummy, 0, NULL, geom0) == NULL,
             "n_atoms<=0 rejected (vib analysis)");
}

int main(void) {
  test_invalid_input();
  test_h2_hessian_and_frequency();

  if (failures == 0) {
    printf("\nAll test_vibrational checks passed.\n");
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
  }

  return failures == 0 ? 0 : 1;
}
