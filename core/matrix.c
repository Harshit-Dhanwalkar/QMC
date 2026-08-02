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
  cmatrix_t *m = malloc(sizeof(cmatrix_t));
  if (!m) {
    return NULL;
  }

  m->data = calloc((size_t)nrows * ncols, sizeof(complex_t));
  if (!m->data) {
    free(m);

    return NULL;
  }

  m->nrows = nrows;
  m->ncols = ncols;

  return m;
}

void cmatrix_free(cmatrix_t *m) {
  if (!m) {
    return;
  }

  free(m->data);
  free(m);
}

cmatrix_t *cmatrix_copy(const cmatrix_t *m) {
  if (!m) {
    return NULL;
  }

  cmatrix_t *c = cmatrix_alloc(m->nrows, m->ncols);
  if (!c) {
    return NULL;
  }

  memcpy(c->data, m->data, (size_t)m->nrows * m->ncols * sizeof(complex_t));

  return c;
}

// Basic operations
cmatrix_t *cmatrix_multiply(const cmatrix_t *a, const cmatrix_t *b) {
  if (!a || !b || a->ncols != b->nrows) {
    return NULL;
  }

  cmatrix_t *c = cmatrix_alloc(a->nrows, b->ncols);
  if (!c) {
    return NULL;
  }

  for (int i = 0; i < a->nrows; i++) {
    for (int j = 0; j < b->ncols; j++) {
      complex_t sum = c_zero();

      for (int k = 0; k < a->ncols; k++) {
        sum = c_add(sum, c_mul(CMAT(a, i, k), CMAT(b, k, j)));
      }

      CMAT(c, i, j) = sum;
    }
  }

  return c;
}

cmatrix_t *cmatrix_transpose(const cmatrix_t *m) {
  if (!m) {
    return NULL;
  }

  cmatrix_t *t = cmatrix_alloc(m->ncols, m->nrows);
  if (!t) {
    return NULL;
  }

  for (int i = 0; i < m->nrows; i++) {
    for (int j = 0; j < m->ncols; j++) {
      CMAT(t, j, i) = CMAT(m, i, j);
    }
  }

  return t;
}

cmatrix_t *cmatrix_adjoint(const cmatrix_t *m) {
  if (!m) {
    return NULL;
  }

  cmatrix_t *a = cmatrix_transpose(m);
  if (!a) {
    return NULL;
  }

  for (int i = 0; i < a->nrows; i++) {
    for (int j = 0; j < a->ncols; j++) {
      CMAT(a, i, j) = c_conj(CMAT(a, i, j));
    }
  }

  return a;
}

void cmatrix_scale(cmatrix_t *m, complex_t s) {
  if (!m) {
    return;
  }

  for (int i = 0; i < m->nrows * m->ncols; i++) {
    m->data[i] = c_mul(m->data[i], s);
  }
}

// Dense matrix-vector product y = A*x.
// Returns a allocated vector, or NULL on dimension mismatch/allocation failure.
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
