/*
Angular momentum operators, Clebsch-Gordan, and spin-1/2,
operators, spherical harmonics, spin
*/

#include "angular.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/special/special.h"
#include "../core/vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Pauli matrices (row-major)
complex_t sigma_x[4] = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}};
complex_t sigma_y[4] = {{0.0, 0.0}, {0.0, -1.0}, {0.0, 1.0}, {0.0, 0.0}};
complex_t sigma_z[4] = {{1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {-1.0, 0.0}};

// Ladder operators
complex_t l_plus_op(int l, int m, int m_prime) {
  if (m_prime != m + 1)
    return c_zero();
  if (m < -l || m >= l)
    return c_zero();
  double val = sqrt((double)(l - m) * (l + m + 1));
  return c_real(val);
}

complex_t l_minus_op(int l, int m, int m_prime) {
  if (m_prime != m - 1)
    return c_zero();
  if (m <= -l || m > l)
    return c_zero();
  double val = sqrt((double)(l + m) * (l - m + 1));
  return c_real(val);
}

// Matrix representations
cmatrix_t *l2_matrix(int l) {
  int dim = 2 * l + 1;
  cmatrix_t *M = cmatrix_alloc(dim, dim);
  if (!M)
    return NULL;
  for (int i = 0; i < dim; i++) {
    CMAT(M, i, i) = c_real(l * (l + 1.0));
  }
  return M;
}

cmatrix_t *lz_matrix(int l) {
  int dim = 2 * l + 1;
  cmatrix_t *M = cmatrix_alloc(dim, dim);
  if (!M)
    return NULL;
  for (int i = 0; i < dim; i++) {
    int m = i - l;
    CMAT(M, i, i) = c_real(m);
  }
  return M;
}

cmatrix_t *lx_matrix(int l) {
  int dim = 2 * l + 1;
  cmatrix_t *M = cmatrix_alloc(dim, dim);
  if (!M)
    return NULL;
  for (int i = 0; i < dim; i++) {
    int m = i - l;
    // <m+1|Lx|m> = 0.5 * sqrt((l-m)(l+m+1))
    if (m < l) {
      double val = 0.5 * sqrt((double)(l - m) * (l + m + 1));
      CMAT(M, i + 1, i) = c_real(val);
      CMAT(M, i, i + 1) = c_real(val);
    }
  }
  return M;
}

cmatrix_t *ly_matrix(int l) {
  int dim = 2 * l + 1;
  cmatrix_t *M = cmatrix_alloc(dim, dim);
  if (!M)
    return NULL;
  for (int i = 0; i < dim; i++) {
    int m = i - l;
    if (m < l) {
      double val = 0.5 * sqrt((double)(l - m) * (l + m + 1));
      CMAT(M, i + 1, i) = c_imag(-val);
      CMAT(M, i, i + 1) = c_imag(val);
    }
  }
  return M;
}

static double lfact(int n) { return lgamma((double)n + 1.0); }

// Clebsch-Gordan coefficient via Racah's formula.
double clebsch_gordan(int j1_2, int m1_2, int j2_2, int m2_2, int J_2,
                      int M_2) {
  if (m1_2 + m2_2 != M_2)
    return 0.0;
  if (J_2 < abs(j1_2 - j2_2) || J_2 > j1_2 + j2_2)
    return 0.0;
  if (abs(m1_2) > j1_2 || abs(m2_2) > j2_2 || abs(M_2) > J_2)
    return 0.0;
  if ((j1_2 + j2_2 + J_2) % 2 != 0)
    return 0.0;

  int j1 = j1_2, m1 = m1_2, j2 = j2_2, m2 = m2_2, J = J_2, M = M_2;

  double log_pref =
      0.5 * (log((double)(J + 1)) + lfact((J + j1 - j2) / 2) +
             lfact((J - j1 + j2) / 2) + lfact((j1 + j2 - J) / 2) -
             lfact((j1 + j2 + J) / 2 + 1) + lfact((J + M) / 2) +
             lfact((J - M) / 2) + lfact((j1 - m1) / 2) + lfact((j1 + m1) / 2) +
             lfact((j2 - m2) / 2) + lfact((j2 + m2) / 2));

  int kmin = 0;
  kmin = kmin > -(J - j2 + m1) / 2 ? kmin : -(J - j2 + m1) / 2;
  kmin = kmin > -(J - j1 - m2) / 2 ? kmin : -(J - j1 - m2) / 2;
  int kmax = (j1 + j2 - J) / 2;
  int t1 = (j1 - m1) / 2, t2 = (j2 + m2) / 2;
  kmax = kmax < t1 ? kmax : t1;
  kmax = kmax < t2 ? kmax : t2;

  double sum = 0.0;
  for (int k = kmin; k <= kmax; k++) {
    double log_term =
        -(lfact(k) + lfact((j1 + j2 - J) / 2 - k) + lfact((j1 - m1) / 2 - k) +
          lfact((j2 + m2) / 2 - k) + lfact((J - j2 + m1) / 2 + k) +
          lfact((J - j1 - m2) / 2 + k));
    double term = exp(log_term);
    if (k % 2 != 0)
      term = -term;
    sum += term;
  }

  return exp(log_pref) * sum;
}

// Spin operations
void spin_sigma_x(cvector_t *spinor) {
  if (!spinor || spinor->n != 2)
    return;
  complex_t a = spinor->data[0];
  spinor->data[0] = spinor->data[1];
  spinor->data[1] = a;
}

void spin_sigma_y(cvector_t *spinor) {
  if (!spinor || spinor->n != 2)
    return;
  complex_t a = spinor->data[0];
  spinor->data[0] = c_mul(c_imag(-1.0), spinor->data[1]);
  spinor->data[1] = c_mul(c_imag(1.0), a);
}

void spin_sigma_z(cvector_t *spinor) {
  if (!spinor || spinor->n != 2)
    return;
  spinor->data[0] = spinor->data[0];
  spinor->data[1] = c_scale(spinor->data[1], -1.0);
}
