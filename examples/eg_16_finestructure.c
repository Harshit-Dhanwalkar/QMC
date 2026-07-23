/*
 * Hydrogen Fine Structure (Spin-Orbit Coupling)
 *
 * Demonstrates spin_orbit_ls_expect_from_coupling (built from
 * couple_states + ladder operators) matching the closed-form
 * j(j+1)-l(l+1)-s(s+1) formula, then known total fine-structure
 * shift, checked against theoretical 2P_{3/2}-2P_{1/2} splitting
 * (~10.97 GHz).
 */

#include "../core/constants.h"
#include "../physics/angular.h"
#include "../physics/fine_structure.h"
#include "../physics/hydrogen.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  printf(" > Hydrogen Fine Structure (Spin-Orbit Coupling)\n\n");

  printf("   <L.S>/\\hbar^2: closed form vs. CG-coupling reconstruction\n");
  printf("   l   j     closed-form   from-coupling   diff\n");
  printf("  ---  ---   -----------   -------------   --------\n");
  for (int l = 0; l <= 2; l++) {
    int J2_list[4];
    int count = couple_allowed_J(2 * l, 1, J2_list);

    for (int idx = 0; idx < count; idx++) {
      int j_2 = J2_list[idx];
      double closed = spin_orbit_ls_expect(l, j_2);
      double coupled = spin_orbit_ls_expect_from_coupling(l, j_2, j_2);
      printf("   %d  %d/2   %11.6f   %13.6f   %8.2e\n", l, j_2, closed, coupled,
             fabs(closed - coupled));
    }
  }
  printf("\n");

  printf(
      "   Hydrogen fine-structure energy shifts (all depend only on n,j):\n");
  printf("   n   j     dE_fs (eV)\n");
  printf("  ---  ---   ------------\n");
  struct {
    int n, j_2;
  } states[] = {{2, 1}, {2, 3}, {3, 1}, {3, 3}, {3, 5}};

  for (int i = 0; i < 5; i++) {
    double dE =
        hydrogen_fine_structure_shift(states[i].n, states[i].j_2, HBAR,
                                      M_ELECTRON, E_CHARGE, EPSILON_0, C_LIGHT);
    printf("   %d   %d/2   %12.6e\n", states[i].n, states[i].j_2,
           dE / E_CHARGE);
  }
  printf("\n");

  double E32 = hydrogen_fine_structure_shift(2, 3, HBAR, M_ELECTRON, E_CHARGE,
                                             EPSILON_0, C_LIGHT);
  double E12 = hydrogen_fine_structure_shift(2, 1, HBAR, M_ELECTRON, E_CHARGE,
                                             EPSILON_0, C_LIGHT);
  double split_GHz = ((E32 - E12) / (2.0 * M_PI * HBAR)) / 1e9;
  printf("   2P_3/2 - 2P_1/2 splitting: %.3f GHz (theoretical: ~10.97 GHz)\n",
         split_GHz);

  return 0;
}
