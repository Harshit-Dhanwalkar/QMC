/*
 * Angular Momentum Coupling
 *
 * Demonstrates couple_states() building coupled |J,M> states from
 * Clebsch-Gordan coefficients: two spin-1/2's -> triplet/singlet, and
 * l=1 (x) s=1/2 -> spin-orbit coupling basis.
 */

#include "../core/complex.h"
#include "../core/vector.h"
#include "../physics/angular.h"
#include <math.h>
#include <stdio.h>

static void print_two_spin_half(int J_2, int M_2, const char *name) {
  cvector_t *v = couple_states(1, 1, J_2, M_2);
  if (!v) {
    printf("   (invalid state)\n");

    return;
  }

  printf("   %-10s (J=%d/2,M=%d/2): |down down>=%6.3f  |down up>=%6.3f  "
         "|up down>=%6.3f  |up up>=%6.3f\n",
         name, J_2, M_2, v->data[0].re, v->data[1].re, v->data[2].re,
         v->data[3].re);

  cvector_free(v);
}

int main(void) {
  printf(" > Angular Momentum Coupling (Clebsch-Gordan)\n\n");

  printf("   Two spin-1/2 particles -> triplet (J=1) + singlet (J=0):\n");
  print_two_spin_half(2, 2, "Triplet M=+1");
  print_two_spin_half(2, 0, "Triplet M=0");
  print_two_spin_half(2, -2, "Triplet M=-1");
  print_two_spin_half(0, 0, "Singlet");
  printf("\n");

  printf("   l=1 (orbital) (x) s=1/2 (spin): allowed J and CG coeff.s\n");
  for (int l = 1; l <= 2; l++) {
    int J2_list[4];
    int count = couple_allowed_J(2 * l, 1, J2_list);
    printf("   l=%d: allowed J (doubled) =", l);
    for (int i = 0; i < count; i++) {
      printf(" %d", J2_list[i]);
    }
    printf("\n");

    for (int idx = 0; idx < count; idx++) {
      int J_2 = J2_list[idx];
      int M_2 = 1; // M=1/2, always valid for any half-integer J here
      cvector_t *v = couple_states(2 * l, 1, J_2, M_2);
      if (!v) {
        continue;
      }

      double norm_sq = 0.0;
      for (int i = 0; i < v->n; i++) {
        norm_sq += v->data[i].re * v->data[i].re;
      }
      printf("     J=%d/2, M=1/2: norm^2=%.6f (should be 1.0)\n", J_2, norm_sq);

      cvector_free(v);
    }
  }

  return 0;
}
