/*
Zeeman effect: Lande g-factor, weak-field energy shift, and
Clebsch-Gordan-coupling cross-check of <Sz>.
*/

#include "zeeman.h"
#include "../core/vector.h"
#include "angular.h"
#include <math.h>

static const double HALF = 0.5;
static const double TWO = 2.0;

double zeeman_lande_g_factor(int l_quantum, int j_doubled) {
  if (j_doubled != 2 * l_quantum + 1 && j_doubled != 2 * l_quantum - 1) {
    return NAN; // not a valid l_quantum (x) spin-1/2 coupling
  }

  if (j_doubled <= 0) {
    return NAN; // j=0: formula divides by j(j+1)=0
  }

  double l = (double)l_quantum;
  double s = HALF;
  double j = j_doubled / 2.0;

  double num = j * (j + 1.0) + s * (s + 1.0) - l * (l + 1.0);
  double den = TWO * j * (j + 1.0);

  return 1.0 + num / den;
}

double zeeman_energy_shift(int l_quantum, int j_doubled, int mj_doubled,
                           double B_field, double mu_Bohr) {
  double g_J = zeeman_lande_g_factor(l_quantum, j_doubled);
  if (isnan(g_J)) {
    return NAN;
  }

  double mj = mj_doubled / TWO;

  return g_J * mu_Bohr * B_field * mj;
}

double zeeman_sz_expect_from_coupling(int l_quantum, int j_doubled,
                                      int mj_doubled) {
  int j1_2 = 2 * l_quantum; // orbital (doubled)
  int j2_2 = 1;             // spin-1/2 (doubled)
  cvector_t *v = couple_states(j1_2, j2_2, j_doubled, mj_doubled);
  if (!v) {
    return NAN;
  }

  int dim1 = j1_2 + 1; // 2l+1
  int dim2 = j2_2 + 1; // 2

  double sz = 0.0;
  for (int i1 = 0; i1 < dim1; i1++) {
    for (int i2 = 0; i2 < dim2; i2++) {
      int m2_2 = -j2_2 + 2 * i2; // doubled spin m: -1 (down) or +1 (up)
      double m2 = m2_2 / TWO;
      double c_here = v->data[i1 * dim2 + i2].re;

      sz += c_here * c_here * m2;
    }
  }

  cvector_free(v);

  return sz;
}
