// TODO: Add docs + docstrings

#include "sparse.h"
#include "complex.h"
#include "linalg/tridiag_eigh.h"
#include "matrix.h"
#include "random.h"
#include "vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Allocation / Free
sparse_matrix_t *sparse_alloc(int nrows, int ncols, int nnz) {
  sparse_matrix_t *A = malloc(sizeof(sparse_matrix_t));
  if (!A) {
    return NULL;
  }

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
  if (!A) {
    return;
  }

  free(A->row_ptr);
  free(A->col_ind);
  free(A->values);
  free(A);
}

// Dense to sparse
sparse_matrix_t *sparse_from_dense(const cmatrix_t *A, double tol) {
  if (!A) {
    return NULL;
  }

  int nrows = A->nrows, ncols = A->ncols;
  int nnz = 0;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      if (c_abs(CMAT(A, i, j)) > tol) {
        nnz++;
      }
    }
  }

  sparse_matrix_t *S = sparse_alloc(nrows, ncols, nnz);
  if (!S) {
    return NULL;
  }

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
  if (!A || !x || !y || A->ncols != x->n || A->nrows != y->n) {
    return;
  }

  for (int i = 0; i < A->nrows; i++) {
    complex_t sum = c_zero();

    for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
      sum = c_add(sum, c_mul(A->values[j], x->data[A->col_ind[j]]));
    }

    y->data[i] = sum;
  }
}

// Lanczos
/*
 * Lanczos iteration for k lowest eigenvalues/eigenvectors of Hermitian sparse
 * matrix A.
 *
 * Three-term Krylov recurrence, starting from a fixed-seed random vector.
 * NOTE: LANCZOS_SEED below: a random start has high-probability nonzero overlap
 * with every eigenvector
 *
 *   work_j      = A v_j - \beta_{j-1} v_{j-1}
 *   \alpha_j = Re(<v_j, work_j>)          (real for Hermitian A)
 *   work_j     -= \alpha_j v_j
 *   \beta_j  = ||wwork_j||,
 *   v_{j+1}  = work_j / \beta_j
 *
 * Uses full reorthogonalization (Gram-Schmidt of work_j against every previous
 * v_i each step).
 *
 * Returns NULL on invalid input (non-square A, k<1, k>n, max_iter<k, tol<=0),
 * on allocation failure, or if Krylov subspace collapses (invariant subspace
 * found) before k directions have been generated.
 */
#define LANCZOS_SEED 0x1A2C205ULL

lanczos_result_t *lanczos_eigs(const sparse_matrix_t *A, int k, int max_iter,
                               double tol) {
  if (!A || A->nrows != A->ncols || k < 1 || max_iter < k || tol <= 0.0) {
    return NULL;
  }

  int rows = A->nrows;
  if (k > rows) {
    return NULL;
  }

  int krylov_dim =
      (max_iter < rows) ? max_iter : rows; // Krylov subspace can't exceed rows

  cvector_t **vecs = calloc((size_t)krylov_dim, sizeof *vecs);
  double *alpha = malloc((size_t)krylov_dim * sizeof *alpha);
  double *beta =
      (krylov_dim > 1) ? malloc((size_t)(krylov_dim - 1) * sizeof *beta) : NULL;
  cvector_t *work = cvector_alloc(rows);
  if (!vecs || !alpha || (krylov_dim > 1 && !beta) || !work) {
    free(vecs);
    free(alpha);
    free(beta);
    cvector_free(work);

    return NULL;
  }

  rng_state_t rng;
  rng_seed(&rng, LANCZOS_SEED);

  vecs[0] = cvector_alloc(rows);
  if (!vecs[0]) {
    free(vecs);
    free(alpha);
    free(beta);
    cvector_free(work);

    return NULL;
  }

  for (int i = 0; i < rows; i++) {
    vecs[0]->data[i] = c_new(rng_gaussian(&rng), rng_gaussian(&rng));
  }
  cvector_normalize(vecs[0]);

  int m_eff = 0;

  for (int j = 0; j < krylov_dim; j++) {
    sparse_mv(A, vecs[j], work);

    if (j > 0) {
      for (int idx = 0; idx < rows; idx++) {
        work->data[idx] = c_sub(work->data[idx],
                                c_scale(vecs[j - 1]->data[idx], beta[j - 1]));
      }
    }

    alpha[j] = cvector_expect(vecs[j], work); // Re(<vecs_j, work>)

    for (int idx = 0; idx < rows; idx++) {
      work->data[idx] =
          c_sub(work->data[idx], c_scale(vecs[j]->data[idx], alpha[j]));
    }

    // Full reorthogonalization against every previous Lanczos vector
    for (int i = 0; i <= j; i++) {
      complex_t proj = cvector_dot(vecs[i], work);

      for (int idx = 0; idx < rows; idx++) {
        work->data[idx] =
            c_sub(work->data[idx], c_mul(proj, vecs[i]->data[idx]));
      }
    }

    double bj = cvector_norm(work);
    m_eff = j + 1;

    if (bj < tol || j == krylov_dim - 1) {
      /* Invariant subspace found (bj ~ 0), or out of allotted steps: stop
       * without generating vecs[j+1]. */
      break;
    }

    beta[j] = bj;
    vecs[j + 1] = cvector_alloc(rows);
    if (!vecs[j + 1]) {
      m_eff = j + 1; // only vecs[0..j] were successfully allocated

      break;
    }

    for (int idx = 0; idx < rows; idx++) {
      vecs[j + 1]->data[idx] = c_scale(work->data[idx], 1.0 / bj);
    }
  }

  lanczos_result_t *res = NULL;

  if (m_eff >= k) {
    eigen_t *T_eig = tridiag_eigh(alpha, beta, m_eff);

    if (T_eig) {
      res = malloc(sizeof *res);
      if (res) {
        res->n = rows;
        res->values = malloc((size_t)k * sizeof *res->values);
        res->vectors = cmatrix_alloc(rows, k);

        if (!res->values || !res->vectors) {
          free(res->values);
          cmatrix_free(res->vectors);
          free(res);

          res = NULL;
        } else {
          for (int i = 0; i < k; i++) {
            res->values[i] = T_eig->eigenvalues[i];

            for (int row = 0; row < rows; row++) {
              CMAT(res->vectors, row, i) = c_zero();
            }

            for (int jb = 0; jb < m_eff; jb++) {
              complex_t coeff = CMAT(T_eig->eigenvectors, jb, i);

              for (int row = 0; row < rows; row++) {
                CMAT(res->vectors, row, i) =
                    c_add(CMAT(res->vectors, row, i),
                          c_mul(coeff, vecs[jb]->data[row]));
              }
            }
          }
        }
      }

      eigen_free(T_eig);
    }
  }

  for (int j = 0; j < krylov_dim; j++) {
    if (vecs[j]) {
      cvector_free(vecs[j]);
    }
  }
  free(vecs);
  free(alpha);
  free(beta);
  cvector_free(work);

  return res;
}

void lanczos_free(lanczos_result_t *res) {
  if (!res) {
    return;
  }

  free(res->values);
  cmatrix_free(res->vectors);
  free(res);
}

lanczos_tridiag_t *lanczos_tridiagonalize(const sparse_matrix_t *A,
                                          const cvector_t *v0, int max_iter,
                                          double tol) {
  if (!A || A->nrows != A->ncols || !v0 || v0->n != A->nrows || max_iter < 1 ||
      tol <= 0.0) {
    return NULL;
  }

  int rows = A->nrows;
  int krylov_dim = (max_iter < rows) ? max_iter : rows;

  double v0_norm = cvector_norm(v0);
  if (fabs(v0_norm - 1.0) > 1e-6) {
    return NULL;
  }

  cvector_t **vecs = calloc((size_t)krylov_dim, sizeof *vecs);
  double *alpha = malloc((size_t)krylov_dim * sizeof *alpha);
  double *beta =
      (krylov_dim > 1) ? malloc((size_t)(krylov_dim - 1) * sizeof *beta) : NULL;
  cvector_t *work = cvector_alloc(rows);
  if (!vecs || !alpha || (krylov_dim > 1 && !beta) || !work) {
    free(vecs);
    free(alpha);
    free(beta);
    cvector_free(work);

    return NULL;
  }

  vecs[0] = cvector_copy(v0);
  if (!vecs[0]) {
    free(vecs);
    free(alpha);
    free(beta);
    cvector_free(work);

    return NULL;
  }

  int m_eff = 0;

  for (int j = 0; j < krylov_dim; j++) {
    sparse_mv(A, vecs[j], work);

    if (j > 0) {
      for (int idx = 0; idx < rows; idx++) {
        work->data[idx] = c_sub(work->data[idx],
                                c_scale(vecs[j - 1]->data[idx], beta[j - 1]));
      }
    }

    alpha[j] = cvector_expect(vecs[j], work); // Re(<vecs_j, work>)

    for (int idx = 0; idx < rows; idx++) {
      work->data[idx] =
          c_sub(work->data[idx], c_scale(vecs[j]->data[idx], alpha[j]));
    }

    // Full reorthogonalization against every previous Lanczos vector
    for (int i = 0; i <= j; i++) {
      complex_t proj = cvector_dot(vecs[i], work);

      for (int idx = 0; idx < rows; idx++) {
        work->data[idx] =
            c_sub(work->data[idx], c_mul(proj, vecs[i]->data[idx]));
      }
    }

    double bj = cvector_norm(work);
    m_eff = j + 1;

    if (bj < tol || j == krylov_dim - 1) {
      break;
    }

    beta[j] = bj;
    vecs[j + 1] = cvector_alloc(rows);
    if (!vecs[j + 1]) {
      m_eff = j + 1;

      break;
    }

    for (int idx = 0; idx < rows; idx++) {
      vecs[j + 1]->data[idx] = c_scale(work->data[idx], 1.0 / bj);
    }
  }

  for (int j = 0; j < krylov_dim; j++) {
    if (vecs[j]) {
      cvector_free(vecs[j]);
    }
  }

  free(vecs);
  cvector_free(work);

  lanczos_tridiag_t *res = malloc(sizeof *res);
  if (!res) {
    free(alpha);
    free(beta);

    return NULL;
  }

  res->m = m_eff;
  res->alpha = realloc(alpha, (size_t)m_eff * sizeof *res->alpha);
  if (!res->alpha) {
    res->alpha = alpha; /* realloc to smaller size should not fail, but if it
                           somehow does, the original block is still valid */
  }

  if (m_eff > 1) {
    res->beta = realloc(beta, (size_t)(m_eff - 1) * sizeof *res->beta);
    if (!res->beta) {
      res->beta = beta;
    }
  } else {
    res->beta = NULL;

    free(beta);
  }

  return res;
}

void lanczos_tridiag_free(lanczos_tridiag_t *t) {
  if (!t) {
    return;
  }

  free(t->alpha);
  free(t->beta);
  free(t);
}
