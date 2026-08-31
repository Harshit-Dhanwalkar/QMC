#include "matrix.h"
#include "complex.h"
#include "linalg/lu.h"
#include "vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Allocation / Free / Copy
cmatrix_t *cmatrix_alloc(int nrows, int ncols) {
  cmatrix_t *matrix = malloc(sizeof(cmatrix_t));
  if (!matrix) {
    return NULL;
  }

  matrix->data = calloc((size_t)nrows * ncols, sizeof(complex_t));
  if (!matrix->data) {
    free(matrix);

    return NULL;
  }

  matrix->nrows = nrows;
  matrix->ncols = ncols;

  return matrix;
}

void cmatrix_free(cmatrix_t *matrix) {
  if (!matrix) {
    return;
  }

  free(matrix->data);
  free(matrix);
}

cmatrix_t *cmatrix_copy(const cmatrix_t *matrix) {
  if (!matrix) {
    return NULL;
  }

  cmatrix_t *copy = cmatrix_alloc(matrix->nrows, matrix->ncols);
  if (!copy) {
    return NULL;
  }

  memcpy(copy->data, matrix->data,
         (size_t)matrix->nrows * matrix->ncols * sizeof(complex_t));

  return copy;
}

// Basic operations
cmatrix_t *cmatrix_multiply(const cmatrix_t *left, const cmatrix_t *right) {
  if (!left || !right || left->ncols != right->nrows) {
    return NULL;
  }

  cmatrix_t *copy = cmatrix_alloc(left->nrows, right->ncols);
  if (!copy) {
    return NULL;
  }

  for (int i = 0; i < left->nrows; i++) {
    for (int j = 0; j < right->ncols; j++) {
      complex_t sum = c_zero();

      for (int k = 0; k < left->ncols; k++) {
        sum = c_add(sum, c_mul(CMAT(left, i, k), CMAT(right, k, j)));
      }

      CMAT(copy, i, j) = sum;
    }
  }

  return copy;
}

cmatrix_t *cmatrix_transpose(const cmatrix_t *matrix) {
  if (!matrix) {
    return NULL;
  }

  cmatrix_t *transpose = cmatrix_alloc(matrix->ncols, matrix->nrows);
  if (!transpose) {
    return NULL;
  }

  for (int i = 0; i < matrix->nrows; i++) {
    for (int j = 0; j < matrix->ncols; j++) {
      CMAT(transpose, j, i) = CMAT(matrix, i, j);
    }
  }

  return transpose;
}

cmatrix_t *cmatrix_adjoint(const cmatrix_t *matrix) {
  if (!matrix) {
    return NULL;
  }

  cmatrix_t *left = cmatrix_transpose(matrix);
  if (!left) {
    return NULL;
  }

  for (int i = 0; i < left->nrows; i++) {
    for (int j = 0; j < left->ncols; j++) {
      CMAT(left, i, j) = c_conj(CMAT(left, i, j));
    }
  }

  return left;
}

void cmatrix_scale(cmatrix_t *matrix, complex_t s) {
  if (!matrix) {
    return;
  }

  for (int i = 0; i < matrix->nrows * matrix->ncols; i++) {
    matrix->data[i] = c_mul(matrix->data[i], s);
  }
}

/* Dense matrix-vector product y = A*x (matches sparse_mv's naming/convention
 * for sparse case).
 *
 * Returns a allocated vector, or NULL on dimension  mismatch/allocation
 * failure.
 */
cvector_t *cmatrix_mv(const cmatrix_t *a, const cvector_t *x) {
  if (!a || !x || a->ncols != x->n) {
    return NULL;
  }

  cvector_t *y = cvector_alloc(a->nrows);
  if (!y) {
    return NULL;
  }

  for (int i = 0; i < a->nrows; i++) {
    complex_t sum = c_zero();

    for (int j = 0; j < a->ncols; j++) {
      sum = c_add(sum, c_mul(CMAT(a, i, j), x->data[j]));
    }

    y->data[i] = sum;
  }

  return y;
}

/* Element-wise matrix sum c = a + b.
 *
 * Returns NULL on dimension mismatch or allocation failure.
 */
cmatrix_t *cmatrix_add(const cmatrix_t *a, const cmatrix_t *b) {
  if (!a || !b || a->nrows != b->nrows || a->ncols != b->ncols) {
    return NULL;
  }

  cmatrix_t *copy = cmatrix_alloc(a->nrows, a->ncols);
  if (!copy) {
    return NULL;
  }

  for (int i = 0; i < a->nrows * a->ncols; i++) {
    copy->data[i] = c_add(a->data[i], b->data[i]);
  }

  return copy;
}

// LU decomposition (wrapper)
void cmatrix_lu_decomp(cmatrix_t *a, int *pivot) {
  if (!a || !pivot) {
    return;
  }

  int *p = lu_decompose(a);
  if (p) {
    memcpy(pivot, p, a->nrows * sizeof(int));

    free(p);
  }
}

// Solve A x = b (uses LU)
cvector_t *cmatrix_solve(cmatrix_t *a, const cvector_t *b) {
  if (!a || !b || a->nrows != a->ncols || a->nrows != b->n) {
    return NULL;
  }

  int *pivot = lu_decompose(a);
  if (!pivot) {
    return NULL;
  }

  cvector_t *x = cvector_alloc(a->nrows);
  if (!x) {
    free(pivot);

    return NULL;
  }

  if (lu_solve(a, pivot, b, x) != 0) {
    cvector_free(x);
    free(pivot);

    return NULL;
  }

  free(pivot);

  return x;
}

// Print
void cmatrix_print(const cmatrix_t *m, const char *label) {
  if (label) {
    printf("%s:\n", label);
  }

  if (!m) {
    printf("NULL matrix\n");

    return;
  }

  for (int i = 0; i < m->nrows; i++) {
    for (int j = 0; j < m->ncols; j++) {
      complex_t z = CMAT(m, i, j);

      printf("(%6.3f %+6.3fi) ", z.re, z.im);
    }

    printf("\n");
  }
}
