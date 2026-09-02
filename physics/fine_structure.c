/*
Hydrogen fine structure: spin-orbit coupling + known closed-form formula.
*/

#include "fine_structure.h"
#include "../core/complex.h"
#include "../core/vector.h"
#include "angular.h"
#include "hydrogen.h"
#include <math.h>
#include <stdlib.h>

static const double HALF = 0.5;
static const double TWO = 2.0;

double spin_orbit_ls_expect(int l, int j_2) {
  double j = j_2 / 2.0;
  double ll = (double)l;
  double s = 0.5;

  return 0.5 * (j * (j + 1.0) - ll * (ll + 1.0) - s * (s + 1.0));
}

double hydrogen_expect_inv_r3(int n, int l, double hbar, double mass,
                              double e_charge, double eps0) {
  if (l < 1 || n < 1) {
    return 0.0; // diverges at l=0
  }

  double a0 = 4.0 * M_PI * eps0 * hbar * hbar / (mass * e_charge * e_charge);
  double nn = (double)n, ll = (double)l;

  return 1.0 / (a0 * a0 * a0 * nn * nn * nn * ll * (ll + 0.5) * (ll + 1.0));
}

double hydrogen_spin_orbit_energy(int n, int l, int j_2, double hbar,
                                  double mass, double e_charge, double eps0,
                                  double c) {
  if (l < 1) {
    return 0.0;
  }

  double inv_r3 = hydrogen_expect_inv_r3(n, l, hbar, mass, e_charge, eps0);
  double LS = spin_orbit_ls_expect(l, j_2); // units of \hbar^2
  double coeff =
      (e_charge * e_charge) / (8.0 * M_PI * eps0 * mass * mass * c * c);

  return coeff * inv_r3 * hbar * hbar * LS;
}

double hydrogen_fine_structure_shift(int n, int j_2, double hbar, double mass,
                                     double e_charge, double eps0, double c) {
  if (n < 1) {
    return 0.0;
  }

  double alpha = (e_charge * e_charge) / (4.0 * M_PI * eps0 * hbar * c);
  double E_n = hydrogen_energy_level(n);
  double j = j_2 / 2.0;
  double nn = (double)n;

  return E_n * (alpha * alpha / (nn * nn)) * (nn / (j + 0.5) - 0.75);
}

double spin_orbit_ls_expect_from_coupling(int l, int j_2, int M_2) {
  int j1_2 = 2 * l, j2_2 = 1; // orbital (doubled), spin-1/2 (doubled)
  cvector_t *v = couple_states(j1_2, j2_2, j_2, M_2);
  if (!v) {
    return NAN;
  }

  int dim1 = j1_2 + 1; // 2l+1
  int dim2 = j2_2 + 1; // 2

  // L.S = Lz*Sz + 1/2*(L+ S- + L- S+), evaluated on CG-coupled state.
  double energy = 0.0;
  for (int i1 = 0; i1 < dim1; i1++) {
    int m_orb = i1 - l; // orbital m
    for (int i2 = 0; i2 < dim2; i2++) {
      int m2_2 = -j2_2 + 2 * i2; // spin m, doubled: -1 (down) or +1 (up)
      double m_spin = m2_2 / TWO;
      double c_here = v->data[i1 * dim2 + i2].re;
      if (c_here == 0.0) {
        continue;
      }

      // Diagonal Lz * Sz
      energy += c_here * c_here * m_orb * m_spin;

      // L+ S-: valid when m_orb can be raised (m_{orb} + 1<=l) and m2_2 can be
      // lowered (only from spin-up, m_{spin} = +1/2, down to -1/2).
      if (m_orb + 1 <= l && m2_2 - 2 >= -j2_2) {
        int i1p = i1 + 1;
        int i2p = i2 - 1;
        double c_prime = v->data[i1p * dim2 + i2p].re;
        if (c_prime != 0.0) {
          double Lplus =
              l_plus_op(l, m_orb, m_orb + 1).re; // <l,m_{orb} + 1|L+|l,m_{orb}>
          double Sminus =
              sqrt((HALF + m_spin) *
                   (HALF - m_spin + 1.0)); // <s,m_{spin} - 1|S-|s,m_{spin}>

          energy += c_here * c_prime * HALF * Lplus * Sminus;
        }
      }

      // L- S+: valid when m_orb can be lowered (m_{orb} - 1>=-l) and m2_2 can
      // be raised (only from spin-down, m_{spin} = -1/2, up to +1/2).
      if (m_orb - 1 >= -l && m2_2 + 2 <= j2_2) {
        int i1m = i1 - 1;
        int i2m = i2 + 1;
        double c_prime = v->data[i1m * dim2 + i2m].re;
        if (c_prime != 0.0) {
          double Lminus =
              l_minus_op(l, m_orb, m_orb - 1).re; // <l,m_{orb} - 1|L-|l,m_{orb}>
          double Splus =
              sqrt((HALF - m_spin) *
                   (HALF + m_spin + 1.0)); // <s,m_{spin} + 1|S+|s,m_{spin}>

          energy += c_here * c_prime * HALF * Lminus * Splus;
        }
      }
    }
  }

  cvector_free(v);

  return energy;
}
