/*
Test: hydrogen fine structure / spin-orbit coupling.

1. spin_orbit_ls_expect_from_coupling (built from couple_states + ladder
   operators) must match the closed-form spin_orbit_ls_expect formula,
   for every allowed j at l=1,2,3, and be independent of M.
2. l=0 special case: <L.S> must be exactly 0 (no orbital momentum to
   couple to spin).
3. hydrogen_fine_structure_shift must reproduce the well-known hydrogen
   2P_{3/2}-2P_{1/2} splitting (~10.97 GHz / ~4.5e-5 eV in the
   literature).
*/

#include "../core/constants.h"
#include "../physics/angular.h"
#include "../physics/fine_structure.h"
#include "../physics/hydrogen.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);
  return err > tol;
}

static int test_ls_expect_matches_coupling(void) {
  int fail = 0;
  double tol = 1e-10;

  for (int l = 0; l <= 3; l++) {
    int J2_list[4];
    int count = couple_allowed_J(2 * l, 1, J2_list);
    for (int idx = 0; idx < count; idx++) {
      int j_2 = J2_list[idx];
      double closed_form = spin_orbit_ls_expect(l, j_2);

      // Check every valid M for this J
      for (int M_2 = -j_2; M_2 <= j_2; M_2 += 2) {
        double from_coupling = spin_orbit_ls_expect_from_coupling(l, j_2, M_2);
        char label[64];
        snprintf(label, sizeof label, "l=%d j_2=%d M_2=%d <L.S>", l, j_2, M_2);
        fail |= check_close(from_coupling, closed_form, tol, label);
      }
    }
  }
  return fail;
}

static int test_l0_zero(void) {
  double val = spin_orbit_ls_expect(0, 1); // l=0, j=1/2 only possibility
  return check_close(val, 0.0, 1e-12, "l=0 <L.S>");
}

static int test_2p_splitting(void) {
  // 2P_{3/2} (j_2=3) minus 2P_{1/2} (j_2=1), n=2.
  double E32 = hydrogen_fine_structure_shift(2, 3, HBAR, M_ELECTRON, E_CHARGE,
                                             EPSILON_0, C_LIGHT);
  double E12 = hydrogen_fine_structure_shift(2, 1, HBAR, M_ELECTRON, E_CHARGE,
                                             EPSILON_0, C_LIGHT);
  double split_J = E32 - E12;
  double split_eV = split_J / E_CHARGE;
  double split_GHz = (split_J / (2.0 * M_PI * HBAR)) / 1e9;

  printf("  2P_3/2 - 2P_1/2 splitting: %.4e eV (%.3f GHz)\n", split_eV,
         split_GHz);
  printf("  literature: ~4.5e-5 eV (~10.97 GHz)\n");

  // Loose tolerance
  int fail = fabs(split_eV - 4.5e-5) / 4.5e-5 > 0.1;
  fail |= fabs(split_GHz - 10.97) / 10.97 > 0.1;
  return fail;
}

int main(void) {
  int failed = 0;

  printf("<L.S> from coupling vs closed form (all l,j,M):\n");
  failed += test_ls_expect_matches_coupling();

  printf("l=0 <L.S> = 0:\n");
  failed += test_l0_zero();

  printf("Hydrogen 2P fine-structure splitting:\n");
  failed += test_2p_splitting();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");
  return 0;
}
