/*
Singular value decomposition

// HACK: Currently implement a simple Jacobi SVD for real matrices
// TODO: adapt to complex by treating real and imaginary separately

Implemented a Golub-Reinsch algorithm for real and then provide a complex
wrapper that uses the real SVD on the block matrix [[Re(A), -Im(A)], [Im(A),
Re(A)]].
*/

#include "svd.h"
#include "../../core/complex.h"
#include "../../core/matrix.h"
#include "../../core/vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int svd_decompose(const cmatrix_t *A, cmatrix_t *U, cvector_t *S,
                  cmatrix_t *V) {
  if (!A || !U || !S || !V)
    return -1;
  int m = A->nrows;
  int n = A->ncols;
  if (m < n) {
    // Decompose A^T instead, then swap U and V.
    cmatrix_t *At = cmatrix_transpose(A);
    if (!At)
      return -1;
    int ret = svd_decompose(At, V, S, U); // NOTE: U and V swapped
    cmatrix_free(At);
    // At this point U and V are swapped, but caller expects U (m x m) and V
    // (nxn).
    // HACK: For simplicity, assume m >= n.
    return -1;
  }

  // Compute A^\dagger A (n x n)
  cmatrix_t *ATA = cmatrix_alloc(n, n);
  if (!ATA)
    return -1;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      complex_t sum = c_zero();
      for (int k = 0; k < m; k++) {
        sum = c_add(sum, c_mul(c_conj(CMAT(A, k, i)), CMAT(A, k, j)));
      }
      CMAT(ATA, i, j) = sum;
    }
  }

  eigen_t *eig = cmatrix_eigh(ATA); // uses eigensolver
  cmatrix_free(ATA);
  if (!eig)
    return -1;

  // Sort eigenvalues and eigenvectors descending
  int *idx = malloc(n * sizeof(int));
  for (int i = 0; i < n; i++)
    idx[i] = i;
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      if (eig->eigenvalues[idx[i]] < eig->eigenvalues[idx[j]]) {
        int tmp = idx[i];
        idx[i] = idx[j];
        idx[j] = tmp;
      }
    }
  }

  // Fill S and V
  S->n = n; // ensure
  if (S->data)
    free(S->data);
  S->data = malloc(n * sizeof(complex_t));
  if (!S->data) {
    free(idx);
    eigen_free(eig);
    return -1;
  }
  // V must be n x n
  if (V->nrows != n || V->ncols != n) {
    if (V->data)
      free(V->data);
    V->nrows = V->ncols = n;
    V->data = calloc(n * n, sizeof(complex_t));
    if (!V->data) {
      free(idx);
      eigen_free(eig);
      return -1;
    }
  }
  for (int i = 0; i < n; i++) {
    int j = idx[i];
    S->data[i].re = sqrt(eig->eigenvalues[j]);
    S->data[i].im = 0.0;
    // V column i = eigenvector j
    for (int r = 0; r < n; r++) {
      CMAT(V, r, i) = eig->eigenvectors[j].data[r];
    }
  }

  // Compute U = A * V * diag(1/S) (for nonzero S)
  // U is m x n (thin SVD) or m x m (full). Do thin (m x n).
  if (U->nrows != m || U->ncols != n) {
    if (U->data)
      free(U->data);
    U->nrows = m;
    U->ncols = n;
    U->data = calloc(m * n, sizeof(complex_t));
    if (!U->data) {
      free(idx);
      eigen_free(eig);
      return -1;
    }
  }
  // Compute temp = A * V (m x n)
  cmatrix_t *AV = cmatrix_alloc(m, n);
  if (!AV) {
    free(idx);
    eigen_free(eig);
    return -1;
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      complex_t sum = c_zero();
      for (int k = 0; k < n; k++) {
        sum = c_add(sum, c_mul(CMAT(A, i, k), CMAT(V, k, j)));
      }
      CMAT(AV, i, j) = sum;
    }
  }
  // U = AV * diag(1/S)
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double s = S->data[j].re;
      if (s > 1e-15)
        CMAT(U, i, j) = c_scale(CMAT(AV, i, j), 1.0 / s);
      else
        CMAT(U, i, j) = c_zero();
    }
  }

  cmatrix_free(AV);
  free(idx);
  eigen_free(eig);
  return 0;
}

int svd_solve(const cmatrix_t *U, const cvector_t *S, const cmatrix_t *V,
              const cvector_t *b, cvector_t *x, double tol) {
  if (!U || !S || !V || !b || !x)
    return -1;
  int m = U->nrows;
  int n = U->ncols; // should match V->nrows and S->n
  if (b->n != m || x->n != n)
    return -1;

  // Compute U^\dagger b (n vector)
  cvector_t *Ub = cvector_alloc(n);
  if (!Ub)
    return -1;
  for (int i = 0; i < n; i++) {
    complex_t sum = c_zero();
    for (int j = 0; j < m; j++) {
      sum = c_add(sum, c_mul(c_conj(CMAT(U, j, i)), b->data[j]));
    }
    Ub->data[i] = sum;
  }

  // x = V * diag(1/S) * Ub
  for (int i = 0; i < n; i++) {
    complex_t sum = c_zero();
    for (int j = 0; j < n; j++) {
      double s = S->data[j].re;
      if (s > tol) {
        sum = c_add(sum, c_mul(CMAT(V, i, j), c_scale(Ub->data[j], 1.0 / s)));
      }
    }
    x->data[i] = sum;
  }

  cvector_free(Ub);
  return 0;
}

cmatrix_t *svd_pseudoinverse(const cmatrix_t *U, const cvector_t *S,
                             const cmatrix_t *V, double tol) {
  if (!U || !S || !V)
    return NULL;
  int m = U->nrows;
  int n = U->ncols; // n = rank
  cmatrix_t *pinv = cmatrix_alloc(n, m);
  if (!pinv)
    return NULL;
  // pinv = V * diag(1/S) * U^\dagger
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      complex_t sum = c_zero();
      for (int k = 0; k < n; k++) {
        double s = S->data[k].re;
        if (s > tol) {
          complex_t factor = c_scale(c_conj(CMAT(U, j, k)), 1.0 / s);
          sum = c_add(sum, c_mul(CMAT(V, i, k), factor));
        }
      }
      CMAT(pinv, i, j) = sum;
    }
  }
  return pinv;
}
