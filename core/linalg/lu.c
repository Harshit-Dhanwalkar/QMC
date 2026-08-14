/*
LU decomposition
*/
#include "lu.h"
#include "../complex.h"
#include "../matrix.h"
#include "../vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int *lu_decompose(cmatrix_t *A) {
  if (!A || A->nrows != A->ncols) {
    return NULL;
  }

  int n = A->nrows;
  int *pivot = malloc(n * sizeof(int));
  if (!pivot) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    pivot[i] = i;
  }

  for (int k = 0; k < n; k++) {
    // Find pivot
    double max_abs = -1.0;
    int max_row = k;

    for (int i = k; i < n; i++) {
      double abs_val = c_abs(CMAT(A, i, k));

      if (abs_val > max_abs) {
        max_abs = abs_val;
        max_row = i;
      }
    }

    if (max_abs < 1e-15) {
      free(pivot);

      return NULL; // singular
    }

    // Swap rows
    if (max_row != k) {
      int tmp = pivot[k];

      pivot[k] = pivot[max_row];
      pivot[max_row] = tmp;

      for (int j = 0; j < n; j++) {
        complex_t t = CMAT(A, k, j);
        CMAT(A, k, j) = CMAT(A, max_row, j);
        CMAT(A, max_row, j) = t;
      }
    }

    // Eliminate
    for (int i = k + 1; i < n; i++) {
      complex_t factor = c_div(CMAT(A, i, k), CMAT(A, k, k));
      CMAT(A, i, k) = factor; // store L factor

      for (int j = k + 1; j < n; j++) {
        CMAT(A, i, j) = c_sub(CMAT(A, i, j), c_mul(factor, CMAT(A, k, j)));
      }
    }
  }

  return pivot;
}

int lu_solve(const cmatrix_t *LU, const int *pivot, const cvector_t *b,
             cvector_t *x) {
  if (!LU || !pivot || !b || !x || LU->nrows != LU->ncols) {
    return -1;
  }

  int n = LU->nrows;
  if (b->n != n || x->n != n) {
    return -1;
  }

  // Forward substitution with permutation: Pb = L y
  cvector_t *y = cvector_alloc(n);
  if (!y) {
    return -1;
  }

  for (int i = 0; i < n; i++) {
    complex_t sum = c_zero();

    for (int j = 0; j < i; j++) {
      sum = c_add(sum, c_mul(CMAT(LU, i, j), y->data[j]));
    }

    y->data[i] = c_sub(b->data[pivot[i]], sum);
  }

  // Back substitution: U x = y
  for (int i = n - 1; i >= 0; i--) {
    complex_t sum = c_zero();

    for (int j = i + 1; j < n; j++) {
      sum = c_add(sum, c_mul(CMAT(LU, i, j), x->data[j]));
    }

    x->data[i] = c_div(c_sub(y->data[i], sum), CMAT(LU, i, i));
  }

  cvector_free(y);

  return 0;
}

complex_t lu_det(const cmatrix_t *LU, const int *pivot) {
  if (!LU || LU->nrows != LU->ncols) {
    return c_zero();
  }

  int n = LU->nrows;
  complex_t det = c_one();
  int sign = 1;

  for (int i = 0; i < n; i++) {
    if (pivot[i] != i) {
      sign = -sign;
    }

    det = c_mul(det, CMAT(LU, i, i));
  }

  det = c_scale(det, (double)sign);

  return det;
}

cmatrix_t *lu_invert(const cmatrix_t *A) {
  if (!A || A->nrows != A->ncols) {
    return NULL;
  }

  int n = A->nrows;
  cmatrix_t *A_copy = cmatrix_copy(A);
  if (!A_copy) {
    return NULL;
  }

  int *pivot = lu_decompose(A_copy);
  if (!pivot) {
    cmatrix_free(A_copy);

    return NULL;
  }

  cmatrix_t *inv = cmatrix_alloc(n, n);
  if (!inv) {
    free(pivot);
    cmatrix_free(A_copy);

    return NULL;
  }

  // Solve for each column of identity
  cvector_t *b = cvector_alloc(n);
  cvector_t *x = cvector_alloc(n);
  if (!b || !x) {
    free(pivot);
    cmatrix_free(A_copy);
    cmatrix_free(inv);

    return NULL;
  }

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      b->data[i] = (i == j) ? c_one() : c_zero();
    }

    if (lu_solve(A_copy, pivot, b, x) != 0) {
      free(pivot);
      cmatrix_free(A_copy);
      cmatrix_free(inv);
      cvector_free(b);
      cvector_free(x);

      return NULL;
    }

    for (int i = 0; i < n; i++) {
      CMAT(inv, i, j) = x->data[i];
    }
  }

  cvector_free(b);
  cvector_free(x);
  free(pivot);
  cmatrix_free(A_copy);

  return inv;
}
