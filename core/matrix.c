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
         (size_t)matrix->nrows * matrix->ncols * sizeof(complex_t)); // NOLINT

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

void cmatrix_scale(cmatrix_t *matrix, complex_t scalar) {
  if (!matrix) {
    return;
  }

  for (int i = 0; i < matrix->nrows * matrix->ncols; i++) {
    matrix->data[i] = c_mul(matrix->data[i], scalar);
  }
}

/*
 * Dense matrix-vector product y = A*x (matches sparse_mv's naming/convention
 * for sparse case).
 *
 * Returns a allocated vector, or NULL on dimension  mismatch/allocation
 * failure.
 */
cvector_t *cmatrix_mv(const cmatrix_t *mat, const cvector_t *vec) {
  if (!mat || !vec || mat->ncols != vec->n) {
    return NULL;
  }

  cvector_t *y = cvector_alloc(mat->nrows);
  if (!y) {
    return NULL;
  }

  for (int i = 0; i < mat->nrows; i++) {
    complex_t sum = c_zero();

    for (int j = 0; j < mat->ncols; j++) {
      sum = c_add(sum, c_mul(CMAT(mat, i, j), vec->data[j]));
    }

    y->data[i] = sum;
  }

  return y;
}

/*
 * Element-wise matrix sum c = a + b.
 *
 * Returns NULL on dimension mismatch or allocation failure.
 */
cmatrix_t *cmatrix_add(const cmatrix_t *mat_a, const cmatrix_t *mat_b) {
  if (!mat_a || !mat_b || mat_a->nrows != mat_b->nrows ||
      mat_a->ncols != mat_b->ncols) {
    return NULL;
  }

  cmatrix_t *copy = cmatrix_alloc(mat_a->nrows, mat_a->ncols);
  if (!copy) {
    return NULL;
  }

  for (int i = 0; i < mat_a->nrows * mat_a->ncols; i++) {
    copy->data[i] = c_add(mat_a->data[i], mat_b->data[i]);
  }

  return copy;
}

// LU decomposition (wrapper)
void cmatrix_lu_decomp(cmatrix_t *mat, int *pivot) {
  if (!mat || !pivot) {
    return;
  }

  int *p = lu_decompose(mat);
  if (p) {
    memcpy(pivot, p, mat->nrows * sizeof(int)); // NOLINT

    free(p);
  }
}

// Solve A x = b (uses LU) [mat vec = vec_b]
cvector_t *cmatrix_solve(cmatrix_t *mat, const cvector_t *vec_b) {
  if (!mat || !vec_b || mat->nrows != mat->ncols || mat->nrows != vec_b->n) {
    return NULL;
  }

  int *pivot = lu_decompose(mat);
  if (!pivot) {
    return NULL;
  }

  cvector_t *vec = cvector_alloc(mat->nrows);
  if (!vec) {
    free(pivot);

    return NULL;
  }

  if (lu_solve(mat, pivot, vec_b, vec) != 0) {
    cvector_free(vec);
    free(pivot);

    return NULL;
  }

  free(pivot);

  return vec;
}

// Set
cvector_t *cmatrix_get_row(const cmatrix_t *matrix, int row) {
  if (!matrix || row < 0 || row >= matrix->nrows) {
    return NULL;
  }

  cvector_t *row_vec = cvector_alloc(matrix->ncols);
  if (!row_vec) {
    return NULL;
  }

  for (int col = 0; col < matrix->ncols; col++) {
    row_vec->data[col] = CMAT(matrix, row, col);
  }

  return row_vec;
}

// Ryser's algorithm for computing the matrix permanent
complex_t cmatrix_permanent(const cmatrix_t *matrix) {
  if (!matrix || matrix->nrows != matrix->ncols) {
    return c_zero();
  }

  int n = matrix->nrows;
  if (n == 0) {
    return c_one();
  }

  complex_t perm = c_zero();
  unsigned long long total_subsets = 1ULL << n;

  for (unsigned long long s = 1; s < total_subsets; s++) {
    int set_bits = __builtin_popcountll(s);
    complex_t row_prod = c_one();

    for (int i = 0; i < n; i++) {
      complex_t col_sum = c_zero();
      for (int j = 0; j < n; j++) {
        if (s & (1ULL << j)) {
          col_sum = c_add(col_sum, CMAT(matrix, i, j));
        }
      }
      row_prod = c_mul(row_prod, col_sum);
    }

    if ((n - set_bits) % 2 == 0) {
      perm = c_add(perm, row_prod);
    } else {
      perm = c_sub(perm, row_prod);
    }
  }

  return perm;
}

// Print
void cmatrix_print(const cmatrix_t *mat, const char *label) {
  if (label) {
    printf("%s:\n", label);
  }

  if (!mat) {
    printf("NULL matrix\n");

    return;
  }

  for (int i = 0; i < mat->nrows; i++) {
    for (int j = 0; j < mat->ncols; j++) {
      complex_t z = CMAT(mat, i, j);

      printf("(%6.3f %+6.3fi) ", z.re, z.im);
    }

    printf("\n");
  }
}
