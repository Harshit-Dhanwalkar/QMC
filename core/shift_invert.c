#include "shift_invert.h"

#include "complex.h"
#include "linalg/lu.h"
#include "matrix.h"
#include "random.h"
#include "sparse.h"
#include "vector.h"

#include <math.h>
#include <stdlib.h>

#define SHIFT_INVERT_SEED 0x5417F1ULL

/* Project `v` (in place) onto orthogonal complement of first `n_cols` columns
 * of `basis` (each assumed already normalized), i.e. classical Gram-Schmidt
 * deflation against previously-found eigenvectors
 */
static void deflate_against_columns(cvector_t *v, const cmatrix_t *basis,
                                    int n_cols) {
  int n = v->n;

  for (int col = 0; col < n_cols; col++) {
    complex_t proj = c_zero();
    for (int idx = 0; idx < n; idx++) {
      proj = c_add(proj, c_mul(c_conj(CMAT(basis, idx, col)), v->data[idx]));
    }
    for (int idx = 0; idx < n; idx++) {
      v->data[idx] = c_sub(v->data[idx], c_mul(proj, CMAT(basis, idx, col)));
    }
  }
}

shift_invert_result_t *shift_invert_eigs(const sparse_matrix_t *A, double sigma,
                                         int k, int max_iter, double tol) {
  if (!A || A->nrows != A->ncols || k < 1 || max_iter < 1 || tol <= 0.0) {
    return NULL;
  }

  int n = A->nrows;
  if (k > n) {
    return NULL;
  }

  // Build dense shifted matrix M = A - \sigma * I from A's CSR storage
  cmatrix_t *M = cmatrix_alloc(n, n);
  if (!M) {
    return NULL;
  }

  for (int row = 0; row < n; row++) {
    for (int idx = A->row_ptr[row]; idx < A->row_ptr[row + 1]; idx++) {
      CMAT(M, row, A->col_ind[idx]) = A->values[idx];
    }
  }

  for (int i = 0; i < n; i++) {
    CMAT(M, i, i) = c_sub(CMAT(M, i, i), c_real(sigma));
  }

  int *pivot = lu_decompose(M);
  if (!pivot) {
    cmatrix_free(M);

    return NULL; // \sigma is exactly an eigenvalue of A
  }

  shift_invert_result_t *res = malloc(sizeof(shift_invert_result_t));
  double *values = malloc((size_t)k * sizeof(double));
  cmatrix_t *vectors = cmatrix_alloc(n, k);
  cvector_t *v = cvector_alloc(n);
  cvector_t *w = cvector_alloc(n);
  cvector_t *Av = cvector_alloc(n);
  complex_t *sort_scratch = malloc((size_t)n * sizeof(complex_t));
  if (!res || !values || !vectors || !v || !w || !Av || !sort_scratch) {
    goto fail;
  }

  rng_state_t rng;
  rng_seed(&rng, SHIFT_INVERT_SEED);

  for (int i = 0; i < k; i++) {
    for (int idx = 0; idx < n; idx++) {
      v->data[idx] = c_new(rng_gaussian(&rng), rng_gaussian(&rng));
    }

    deflate_against_columns(v, vectors, i);
    cvector_normalize(v);

    double lambda = sigma;
    int converged = 0;

    for (int it = 0; it < max_iter && !converged; it++) {
      if (lu_solve(M, pivot, v, w) != 0) {
        goto fail;
      }

      deflate_against_columns(w, vectors, i);
      cvector_normalize(w);

      sparse_mv(A, w, Av);
      lambda = cvector_dot(w, Av).re; // Rayleigh quotient wrt original A

      complex_t overlap = cvector_dot(v, w); // <v_prev|w>
      if (fabs(1.0 - c_abs(overlap)) < tol) {
        converged = 1;
      }

      for (int idx = 0; idx < n; idx++) {
        v->data[idx] = w->data[idx];
      }
    }

    values[i] = lambda;
    for (int idx = 0; idx < n; idx++) {
      CMAT(vectors, idx, i) = v->data[idx];
    }
  }

  /* Insertion sort by eigenvalue ascending, permuting matching eigenvector
   * column alongside each key (k is always small - an interior-eigenvalue
   * solver for dense-affordable systems) */
  for (int i = 1; i < k; i++) {
    double key = values[i];

    for (int idx = 0; idx < n; idx++) {
      sort_scratch[idx] = CMAT(vectors, idx, i);
    }

    int j = i - 1;
    while (j >= 0 && values[j] > key) {
      values[j + 1] = values[j];

      for (int idx = 0; idx < n; idx++) {
        CMAT(vectors, idx, j + 1) = CMAT(vectors, idx, j);
      }
      j--;
    }

    values[j + 1] = key;
    for (int idx = 0; idx < n; idx++) {
      CMAT(vectors, idx, j + 1) = sort_scratch[idx];
    }
  }

  res->n = k;
  res->values = values;
  res->vectors = vectors;

  cvector_free(v);
  cvector_free(w);
  cvector_free(Av);
  free(sort_scratch);
  free(pivot);
  cmatrix_free(M);

  return res;

fail:
  free(res);
  free(values);
  cmatrix_free(vectors);
  cvector_free(v);
  cvector_free(w);
  cvector_free(Av);
  free(sort_scratch);
  free(pivot);
  cmatrix_free(M);

  return NULL;
}

void shift_invert_result_free(shift_invert_result_t *res) {
  if (!res) {
    return;
  }

  free(res->values);
  cmatrix_free(res->vectors);
  free(res);
}
