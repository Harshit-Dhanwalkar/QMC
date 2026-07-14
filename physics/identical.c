/*
Identical particles: Slater determinants and symmetrization.
*/

#include "identical.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

cmatrix_t *slater_matrix(cvector_t **orbitals, int N, const int *indices) {
  if (!orbitals || !indices || N < 1)
    return NULL;

  cmatrix_t *M = cmatrix_alloc(N, N);
  if (!M)
    return NULL;

  for (int i = 0; i < N; i++) {
    if (!orbitals[i] || indices[i] < 0 || indices[i] >= orbitals[i]->n) {
      cmatrix_free(M);
      return NULL;
    }
  }

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      CMAT(M, i, j) = orbitals[i]->data[indices[j]];
  return M;
}

// Recursive Laplace (cofactor) expansion along first row.
static complex_t cofactor_expand(cmatrix_t *A, int use_sign) {
  int n = A->nrows;
  if (n == 1)
    return CMAT(A, 0, 0);

  complex_t sum = c_zero();
  double sign = 1.0;
  for (int j = 0; j < n; j++) {
    cmatrix_t *minor = cmatrix_alloc(n - 1, n - 1);
    for (int r = 1; r < n; r++) {
      int mc = 0;
      for (int c = 0; c < n; c++) {
        if (c == j)
          continue;
        CMAT(minor, r - 1, mc) = CMAT(A, r, c);
        mc++;
      }
    }
    complex_t sub = cofactor_expand(minor, use_sign);
    complex_t term = c_mul(CMAT(A, 0, j), sub);
    if (use_sign)
      term = c_scale(term, sign);

    sum = c_add(sum, term);
    cmatrix_free(minor);
    sign = -sign;
  }
  return sum;
}

static double factorial(int n) {
  double r = 1.0;
  for (int i = 2; i <= n; i++)
    r *= i;
  return r;
}

complex_t slater_determinant_value(cvector_t **orbitals, int N,
                                   const int *indices) {
  cmatrix_t *M = slater_matrix(orbitals, N, indices);
  if (!M)
    return c_zero();
  complex_t det = cofactor_expand(M, 1);
  cmatrix_free(M);
  double inv_sqrt_nfact = 1.0 / sqrt(factorial(N));
  return c_scale(det, inv_sqrt_nfact);
}

complex_t bosonic_permanent_value(cvector_t **orbitals, int N,
                                  const int *indices) {
  cmatrix_t *M = slater_matrix(orbitals, N, indices);
  if (!M)
    return c_zero();
  complex_t perm = cofactor_expand(M, 0);
  cmatrix_free(M);
  double inv_sqrt_nfact = 1.0 / sqrt(factorial(N));
  return c_scale(perm, inv_sqrt_nfact);
}
