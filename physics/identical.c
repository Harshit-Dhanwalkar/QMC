/*
Identical particles: Slater determinants and symmetrization.
*/

#include "identical.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include <stdlib.h>
#include <string.h>

cmatrix_t *slater_determinant(cvector_t **orbitals, int N) {
  if (!orbitals || N < 1)
    return NULL;
  // Each orbital is vector of length M (grid points).
  // For an M x N matrix.
  int M = orbitals[0]->n;
  cmatrix_t *mat = cmatrix_alloc(M, N);
  if (!mat)
    return NULL;
  for (int i = 0; i < N; i++) {
    if (orbitals[i]->n != M) {
      cmatrix_free(mat);
      return NULL;
    }
    for (int j = 0; j < M; j++) {
      CMAT(mat, j, i) = orbitals[i]->data[j];
    }
  }
  return mat;
}

cmatrix_t *bosonic_symmetrize(cvector_t **orbitals, int N) {
  // For bosons, permanent is same as the matrix of orbitals,
  // but with normalization factor 1/sqrt(N!).
  // HACK: return the same matrix (normalize yourself).
  return slater_determinant(orbitals, N);
}
