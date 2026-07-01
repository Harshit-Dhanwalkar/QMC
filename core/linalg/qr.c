/*
General QR decomposition (Householder reflections)
*/
#include "qr.h"
#include "../../core/complex.h"
#include "../../core/matrix.h"
#include "../../core/vector.h"
#include "../complex.h"
#include "../matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void householder(cmatrix_t *A, int k, cvector_t *v) {
  // v is length m = n - k
  int n = A->ncols;
  int m = n - k;

  // Extract column k from row k downward
  complex_t *x = malloc(m * sizeof(complex_t));
  for (int i = 0; i < m; i++) {
    x[i] = CMAT(A, k + i, k);
  }

  double norm = 0.0;
  for (int i = 0; i < m; i++)
    norm += c_abs2(x[i]);
  norm = sqrt(norm);
  if (norm < 1e-15) {
    free(x);
    v->n = 0;
    return;
  }

  complex_t alpha = c_neg(x[0]);
  if (c_abs(alpha) > 1e-15) {
    alpha = c_scale(alpha, norm / c_abs(alpha));
  } else {
    alpha = c_real(norm);
  }

  x[0] = c_sub(x[0], alpha);
  double v_norm = 0.0;
  for (int i = 0; i < m; i++)
    v_norm += c_abs2(x[i]);
  v_norm = sqrt(v_norm);
  if (v_norm < 1e-15) {
    free(x);
    v->n = 0;
    return;
  }

  for (int i = 0; i < m; i++)
    x[i] = c_scale(x[i], 1.0 / v_norm);
  // Store v in A (lower part)
  for (int i = 0; i < m; i++) {
    CMAT(A, k + i, k) = x[i];
  }

  // Apply reflection to A (only columns k..n-1)
  for (int j = k; j < n; j++) {
    complex_t sum = c_zero();
    for (int i = 0; i < m; i++) {
      sum = c_add(sum, c_mul(c_conj(x[i]), CMAT(A, k + i, j)));
    }
    complex_t factor = c_scale(sum, -2.0);
    for (int i = 0; i < m; i++) {
      CMAT(A, k + i, j) = c_add(CMAT(A, k + i, j), c_mul(factor, x[i]));
    }
  }

  free(x);
  // Store v in output vector? Not needed; we keep in A.
}

int qr_decompose(cmatrix_t *A, cmatrix_t *Q) {
  if (!A || A->nrows != A->ncols)
    return -1; // HACK: must be square for simplicity; can extend to rectangular

  int n = A->nrows;
  // Initialize Q to identity
  if (Q) {
    if (Q->nrows != n || Q->ncols != n)
      return -1;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        CMAT(Q, i, j) = (i == j) ? c_one() : c_zero();
  }

  // cvector_t *v = cvector_alloc(n);
  for (int k = 0; k < n - 1; k++) {
    // Compute Householder vector in A (modifies A, stores v in column k below
    // diag)
    // HACK: Just call householder which modifies A and also can store v, but
    // need to apply to Q if Q is given. Manually apply to both A and Q.
    int m = n - k;
    // Extract column k from row k downward
    complex_t *x = malloc(m * sizeof(complex_t));
    if (!x)
      return -1;
    for (int i = 0; i < m; i++)
      x[i] = CMAT(A, k + i, k);

    double norm = 0.0;
    for (int i = 0; i < m; i++)
      norm += c_abs2(x[i]);
    norm = sqrt(norm);
    if (norm < 1e-15) {
      free(x);
      continue;
    }

    // Householder vector: v = x - \alpha * e1, \alpha = -sign(x[0]) * norm
    complex_t alpha = c_neg(x[0]);
    if (c_abs(alpha) > 1e-15)
      alpha = c_scale(alpha, norm / c_abs(alpha));
    else
      alpha = c_real(norm);
    x[0] = c_sub(x[0], alpha);

    double v_norm = 0.0;
    for (int i = 0; i < m; i++)
      v_norm += c_abs2(x[i]);
    v_norm = sqrt(v_norm);
    if (v_norm < 1e-15) {
      free(x);
      continue;
    }

    for (int i = 0; i < m; i++)
      x[i] = c_scale(x[i], 1.0 / v_norm);

    // Store Householder vector v in A (lower part)
    for (int i = 0; i < m; i++)
      CMAT(A, k + i, k) = x[i];

    // Apply reflection to A on the right: A = H A, (columns k..n-1)
    for (int j = k; j < n; j++) {
      complex_t sum = c_zero();
      for (int i = 0; i < m; i++) {
        sum = c_add(sum, c_mul(c_conj(x[i]), CMAT(A, k + i, j)));
      }
      complex_t factor = c_scale(sum, -2.0);
      for (int i = 0; i < m; i++) {
        CMAT(A, k + i, j) = c_add(CMAT(A, k + i, j), c_mul(factor, x[i]));
      }
    }

    // Apply reflection to Q (right multiplication: Q = Q^\dagger)
    if (Q) {
      for (int i = 0; i < n; i++) {
        complex_t sum = c_zero();
        for (int j = 0; j < m; j++) {
          sum = c_add(sum, c_mul(CMAT(Q, i, k + j), c_conj(x[j])));
        }

        complex_t factor = c_scale(sum, -2.0);
        for (int j = 0; j < m; j++) {
          CMAT(Q, i, k + j) = c_add(CMAT(Q, i, k + j), c_mul(factor, x[j]));
        }
      }
    }

    free(x);
  }
  // Now A contains R in upper triangle and Householder vectors in lower.
  return 0;
}

int qr_solve(const cmatrix_t *A, const cmatrix_t *Q, const cvector_t *b,
             cvector_t *x) {
  if (!A || !Q || !b || !x)
    return -1;

  int n = A->nrows;
  if (n != A->ncols || n != Q->nrows || n != Q->ncols || b->n != n || x->n != n)
    return -1;

  // Compute Q^\dagger b
  cvector_t *b_tilde = cvector_alloc(n);
  if (!b_tilde)
    return -1;
  for (int i = 0; i < n; i++) {
    complex_t sum = c_zero();
    for (int j = 0; j < n; j++) {
      sum = c_add(sum, c_mul(c_conj(CMAT(Q, j, i)),
                             b->data[j])); // Q^\dagger = conj(Q^T)
    }
    b_tilde->data[i] = sum;
  }

  // Back substitution: R x = b_tilde, (R is upper triangular in A)
  for (int i = n - 1; i >= 0; i--) {
    complex_t sum = c_zero();
    for (int j = i + 1; j < n; j++) {
      sum = c_add(sum, c_mul(CMAT(A, i, j), x->data[j]));
    }
    x->data[i] = c_div(c_sub(b_tilde->data[i], sum), CMAT(A, i, i));
  }

  cvector_free(b_tilde);
  return 0;
}

cmatrix_t *qr_invert(const cmatrix_t *A) {
  if (!A || A->nrows != A->ncols)
    return NULL;

  int n = A->nrows;
  cmatrix_t *A_copy = cmatrix_copy(A);
  if (!A_copy)
    return NULL;

  cmatrix_t *Q = cmatrix_alloc(n, n);
  if (!Q) {
    cmatrix_free(A_copy);
    return NULL;
  }

  if (qr_decompose(A_copy, Q) != 0) {
    cmatrix_free(A_copy);
    cmatrix_free(Q);
    return NULL;
  }

  cmatrix_t *inv = cmatrix_alloc(n, n);
  if (!inv) {
    cmatrix_free(A_copy);
    cmatrix_free(Q);
    return NULL;
  }

  cvector_t *b = cvector_alloc(n);
  cvector_t *x = cvector_alloc(n);
  if (!b || !x) {
    cmatrix_free(A_copy);
    cmatrix_free(Q);
    cmatrix_free(inv);
    return NULL;
  }

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++)
      b->data[i] = (i == j) ? c_one() : c_zero();
    if (qr_solve(A_copy, Q, b, x) != 0) {
      cmatrix_free(A_copy);
      cmatrix_free(Q);
      cmatrix_free(inv);
      cvector_free(b);
      cvector_free(x);
      return NULL;
    }
    for (int i = 0; i < n; i++)
      CMAT(inv, i, j) = x->data[i];
  }

  cvector_free(b);
  cvector_free(x);
  cmatrix_free(A_copy);
  cmatrix_free(Q);
  return inv;
}
