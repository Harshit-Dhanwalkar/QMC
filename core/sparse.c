// TODO: add what it does here

#include "sparse.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Allocation / Free
sparse_matrix_t *sparse_alloc(int nrows, int ncols, int nnz) {
  sparse_matrix_t *A = malloc(sizeof(sparse_matrix_t));
  if (!A)
    return NULL;
  A->nrows = nrows;
  A->ncols = ncols;
  A->nnz = nnz;
  A->row_ptr = malloc((nrows + 1) * sizeof(int));
  A->col_ind = malloc(nnz * sizeof(int));
  A->values = malloc(nnz * sizeof(complex_t));
  if (!A->row_ptr || !A->col_ind || !A->values) {
    sparse_free(A);
    return NULL;
  }
  return A;
}

void sparse_free(sparse_matrix_t *A) {
  if (!A)
    return;
  free(A->row_ptr);
  free(A->col_ind);
  free(A->values);
  free(A);
}

// Dense to sparse
sparse_matrix_t *sparse_from_dense(const cmatrix_t *A, double tol) {
  if (!A)
    return NULL;
  int nrows = A->nrows, ncols = A->ncols;
  int nnz = 0;
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < ncols; j++)
      if (c_abs(CMAT(A, i, j)) > tol)
        nnz++;

  sparse_matrix_t *S = sparse_alloc(nrows, ncols, nnz);
  if (!S)
    return NULL;

  int idx = 0;
  S->row_ptr[0] = 0;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      if (c_abs(CMAT(A, i, j)) > tol) {
        S->col_ind[idx] = j;
        S->values[idx] = CMAT(A, i, j);
        idx++;
      }
    }
    S->row_ptr[i + 1] = idx;
  }
  return S;
}

// Matrix-vector multiply
void sparse_mv(const sparse_matrix_t *A, const cvector_t *x, cvector_t *y) {
  if (!A || !x || !y || A->ncols != x->n || A->nrows != y->n)
    return;
  for (int i = 0; i < A->nrows; i++) {
    complex_t sum = c_zero();
    for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
      sum = c_add(sum, c_mul(A->values[j], x->data[A->col_ind[j]]));
    }
    y->data[i] = sum;
  }
}

// Lanczos
lanczos_result_t *lanczos_eigs(sparse_matrix_t *A, int k, int max_iter,
                               double tol) {
  // TODO: Placeholder: currenlty returns NULL; implement when needed.
  (void)A;
  (void)k;
  (void)max_iter;
  (void)tol;
  return NULL;
}

void lanczos_free(lanczos_result_t *res) {
  if (!res)
    return;
  free(res->values);
  cmatrix_free(res->vectors);
  free(res);
}
