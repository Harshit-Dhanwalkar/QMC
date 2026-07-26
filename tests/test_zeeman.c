/*
Test: zeeman.c (Lande g-factor, weak-field Zeeman shift, and a
Clebsch-Gordan cross-check of <Sz>).

1. Lande g-factors reproduce known textbook values for S1/2, P1/2, P3/2,
   D3/2, D5/2.
2. zeeman_sz_expect_from_coupling (built from couple_states + CG probabilities)
   must match closed-form identity <Sz>=(g_J-1)*mj, for every allowed mj at l=0..3.
3. zeeman_energy_shift must equal g_J * mu_B * B * mj (formula check).
*/

#include "../physics/angular.h"
#include "../physics/zeeman.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

static int test_known_g_factors(void) {
  int fail = 0;
  fail |= check_close(zeeman_lande_g_factor(0, 1), 2.0, 1e-10, "S1/2 g_J");
  fail |=
      check_close(zeeman_lande_g_factor(1, 1), 2.0 / 3.0, 1e-10, "P1/2 g_J");
  fail |=
      check_close(zeeman_lande_g_factor(1, 3), 4.0 / 3.0, 1e-10, "P3/2 g_J");
  fail |=
      check_close(zeeman_lande_g_factor(2, 3), 4.0 / 5.0, 1e-10, "D3/2 g_J");
  fail |=
      check_close(zeeman_lande_g_factor(2, 5), 6.0 / 5.0, 1e-10, "D5/2 g_J");

  return fail;
}

static int test_sz_matches_coupling(void) {
  int fail = 0;
  double tol = 1e-10;

  for (int l = 0; l <= 3; l++) {
    int J2_list[4];
    int count = couple_allowed_J(2 * l, 1, J2_list);

    for (int idx = 0; idx < count; idx++) {
      int j_2 = J2_list[idx];
      double g = zeeman_lande_g_factor(l, j_2);
      if (isnan(g))
        continue;

      for (int mj_2 = -j_2; mj_2 <= j_2; mj_2 += 2) {
        double closed = (g - 1.0) * (mj_2 / 2.0);
        double from_coupling = zeeman_sz_expect_from_coupling(l, j_2, mj_2);

        char label[64];
        snprintf(label, sizeof label, "l=%d j_2=%d mj_2=%d <Sz>", l, j_2, mj_2);
        fail |= check_close(from_coupling, closed, tol, label);
      }
    }
  }

  return fail;
}

static int test_energy_shift_formula(void) {
  int l = 1, j_2 = 3, mj_2 = 3; // P3/2, mj=3/2
  double B = 2.0, mu_B = 1.5;
  double g = zeeman_lande_g_factor(l, j_2);
  double expected = g * mu_B * B * (mj_2 / 2.0);
  double got = zeeman_energy_shift(l, j_2, mj_2, B, mu_B);

  return check_close(got, expected, 1e-12, "dE = g_J * mu_B * B * mj");
}

int main(void) {
  int failed = 0;

  printf("Known Lande g-factors (S1/2,P1/2,P3/2,D3/2,D5/2):\n");
  failed += test_known_g_factors();

  printf("<Sz> from coupling vs closed form (l=0..3, all j, mj):\n");
  failed += test_sz_matches_coupling();

  printf("Energy shift formula:\n");
  failed += test_energy_shift_formula();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
