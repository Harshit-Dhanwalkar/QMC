/*
Zeeman effect: Lande g-factor, weak-field energy shift, and
Clebsch-Gordan-coupling cross-check of <Sz>.
*/

#include "zeeman.h"
#include "../core/vector.h"
#include "angular.h"
#include <math.h>

double zeeman_lande_g_factor(int l, int j_2) {
  if (j_2 != 2 * l + 1 && j_2 != 2 * l - 1) {
    return NAN; // not a valid l (x) spin-1/2 coupling
  }

  if (j_2 <= 0) {
    return NAN; // j=0: formula divides by j(j+1)=0
  }

  double ll = (double)l;
  double s = 0.5;
  double j = j_2 / 2.0;

  double num = j * (j + 1.0) + s * (s + 1.0) - ll * (ll + 1.0);
  double den = 2.0 * j * (j + 1.0);

  return 1.0 + num / den;
}

double zeeman_energy_shift(int l, int j_2, int mj_2, double B, double mu_B) {
  double g_J = zeeman_lande_g_factor(l, j_2);
  if (isnan(g_J)) {
    return NAN;
  }

  double mj = mj_2 / 2.0;

  return g_J * mu_B * B * mj;
}

double zeeman_sz_expect_from_coupling(int l, int j_2, int mj_2) {
  int j1_2 = 2 * l, j2_2 = 1; // orbital (doubled), spin-1/2 (doubled)
  cvector_t *v = couple_states(j1_2, j2_2, j_2, mj_2);
  if (!v) {
    return NAN;
  }

  int dim1 = j1_2 + 1; // 2l+1
  int dim2 = j2_2 + 1; // 2

  double sz = 0.0;
  for (int i1 = 0; i1 < dim1; i1++) {
    for (int i2 = 0; i2 < dim2; i2++) {
      int m2_2 = -j2_2 + 2 * i2; // doubled spin m: -1 (down) or +1 (up)
      double m2 = m2_2 / 2.0;
      double c_here = v->data[i1 * dim2 + i2].re;

      sz += c_here * c_here * m2;
    }
  }

  cvector_free(v);

  return sz;
}
