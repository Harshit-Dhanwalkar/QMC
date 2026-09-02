/*
Complex-Hermitian eigensolver via real-embedding
*/

#include "complex_eigh.h"
#include "../complex.h"
#include "../matrix.h"
#include "eigen_generic.h"
#include <math.h>
#include <stdlib.h>

#ifdef USE_LAPACK
#include <lapacke.h>

// Native complex Hermitian solve via LAPACK's zheev directly on n x n matrix
static eigen_t *cmatrix_eigh_complex_lapack(cmatrix_t *H) {
  int n = H->nrows;

  eigen_t *result = malloc(sizeof(eigen_t));
  if (!result) {
    return NULL;
  }

  result->n = n;
  result->eigenvalues = malloc((size_t)n * sizeof(double));
  result->eigenvectors = cmatrix_alloc(n, n);
  if (!result->eigenvalues || !result->eigenvectors) {
    free(result->eigenvalues);

    if (result->eigenvectors) {
      cmatrix_free(result->eigenvectors);
    }

    free(result);

    return NULL;
  }

  for (int i = 0; i < n * n; i++) {
    result->eigenvectors->data[i] = H->data[i];
  }

  int info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', n,
                           (lapack_complex_double *)result->eigenvectors->data,
                           n, result->eigenvalues);
  if (info != 0) {
    cmatrix_free(result->eigenvectors);

    free(result->eigenvalues);
    free(result);

    return NULL;
  }

  return result;
}
#endif

eigen_t *cmatrix_eigh_complex(cmatrix_t *H) {
  if (!H || H->nrows != H->ncols) {
    return NULL;
  }

#ifdef USE_LAPACK
  return cmatrix_eigh_complex_lapack(H);
#else

  int n = H->nrows;
  int m2 = 2 * n;

  cmatrix_t *M = cmatrix_alloc(m2, m2);
  if (!M) {
    return NULL;
  }

  for (int i = 0; i < m2; i++) {
    for (int j = 0; j < m2; j++) {
      CMAT(M, i, j) = c_zero();
    }
  }

  // M = [[A,-B],[B,A]], H = A + iB
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      complex_t h = CMAT(H, i, j);

      CMAT(M, i, j) = c_real(h.re);
      CMAT(M, i, j + n) = c_real(-h.im);
      CMAT(M, i + n, j) = c_real(h.im);
      CMAT(M, i + n, j + n) = c_real(h.re);
    }
  }

  eigen_t *eig2n = cmatrix_eigh_generic(M);
  cmatrix_free(M);

  if (!eig2n) {
    return NULL;
  }

  if (eig2n->n != m2) {
    eigen_free(eig2n);

    return NULL;
  }

  eigen_t *result = malloc(sizeof *result);

  if (!result) {
    eigen_free(eig2n);

    return NULL;
  }

  result->n = n;
  result->eigenvalues = malloc(n * sizeof(double));
  result->eigenvectors = cmatrix_alloc(n, n);
  if (!result->eigenvalues || !result->eigenvectors) {
    free(result->eigenvalues);

    if (result->eigenvectors) {
      cmatrix_free(result->eigenvectors);
    }

    free(result);
    eigen_free(eig2n);

    return NULL;
  }

  // for (int k = 0; k < n; k++) {
  //   int idx = 2 * k;
  //   double lambda = eig2n->eigenvalues[idx];
  // if (idx + 1 < m2) {
  //   lambda = 0.5 * (lambda + eig2n->eigenvalues[idx + 1]);
  // }
  //
  // result->eigenvalues[k] = lambda;
  //
  // for (int i = 0; i < n; i++) {
  //   double x = CMAT(eig2n->eigenvectors, i, idx).re;
  //   double y = CMAT(eig2n->eigenvectors, i + n, idx).re;
  //   CMAT(result->eigenvectors, i, k) = c_add(c_real(x), c_imag(y));
  // }

  int *order = malloc((size_t)m2 * sizeof(int));

  if (!order) {
    cmatrix_free(result->eigenvectors);
    free(result->eigenvalues);
    free(result);
    eigen_free(eig2n);

    return NULL;
  }

  for (int t = 0; t < m2; t++) {
    order[t] = t;
  }

  for (int a = 1; a < m2; a++) {
    int key = order[a];
    double keyval = eig2n->eigenvalues[key];
    int b = a - 1;

    while (b >= 0 && eig2n->eigenvalues[order[b]] > keyval) {
      order[b + 1] = order[b];
      b--;
    }

    order[b + 1] = key;
  }

  const double degeneracy_tol_rel = 1e-12;
  double spectrum_scale = 0.0;
  for (int t = 0; t < m2; t++) {
    double a = fabs(eig2n->eigenvalues[t]);
    if (a > spectrum_scale) {
      spectrum_scale = a;
    }
  }

  const double degeneracy_tol_floor =
      100.0 * 2.220446049250313e-16 * spectrum_scale;
  const double pairing_tol_factor = 100.0;

  int n_pairs = m2 / 2;
  double *pair_lambda = malloc((size_t)n_pairs * sizeof *pair_lambda);
  int *pair_col0 = malloc((size_t)n_pairs * sizeof *pair_col0);
  int *pair_col1 = malloc((size_t)n_pairs * sizeof *pair_col1);
  if (!pair_lambda || !pair_col0 || !pair_col1) {
    free(pair_lambda);
    free(pair_col0);
    free(pair_col1);
    free(order);
    cmatrix_free(result->eigenvectors);
    free(result->eigenvalues);
    free(result);
    eigen_free(eig2n);

    return NULL;
  }
  for (int p = 0; p < n_pairs; p++) {
    int i0 = order[2 * p], i1 = order[2 * p + 1];
    double v0 = eig2n->eigenvalues[i0], v1 = eig2n->eigenvalues[i1];
    pair_lambda[p] = 0.5 * (v0 + v1);
    pair_col0[p] = i0;
    pair_col1[p] = i1;
    (void)pairing_tol_factor;
  }

  int k_out = 0;
  int pidx = 0;

  while (pidx < n_pairs && k_out < n) {
    int group_start = pidx;
    int group_end = pidx + 1;
    double ref = pair_lambda[group_start];

    const int max_group_pairs = 16;
    while (group_end < n_pairs && group_end - group_start < max_group_pairs) {
      double hi = pair_lambda[group_end];
      double local_scale = fabs(ref) > fabs(hi) ? fabs(ref) : fabs(hi);
      double gap_tol = degeneracy_tol_rel * local_scale + degeneracy_tol_floor;

      if (fabs(hi - ref) >= gap_tol) {
        break;
      }

      group_end++;
    }

    int m_needed = group_end - group_start;
    double lambda_sum = 0.0;

    for (int p = group_start; p < group_end; p++) {
      lambda_sum += pair_lambda[p];
    }

    double lambda = lambda_sum / m_needed;
    int cluster_size = 2 * m_needed;

    // Complex candidate vectors z_t = x_t + i*y_t for every real eigenvector in
    // this cluster
    complex_t *z = malloc((size_t)cluster_size * n * sizeof(complex_t));
    if (!z) {
      free(pair_lambda);
      free(pair_col0);
      free(pair_col1);
      free(order);
      cmatrix_free(result->eigenvectors);
      free(result->eigenvalues);
      free(result);
      eigen_free(eig2n);

      return NULL;
    }

    for (int p = group_start; p < group_end; p++) {
      int t = 2 * (p - group_start);
      const int cols[2] = {pair_col0[p], pair_col1[p]};

      for (int c = 0; c < 2; c++) {
        int col = cols[c];

        for (int i = 0; i < n; i++) {
          double x = CMAT(eig2n->eigenvectors, i, col).re;
          double y = CMAT(eig2n->eigenvectors, i + n, col).re;

          z[(t + c) * n + i] = c_add(c_real(x), c_imag(y));
        }
      }
    }

    /* Modified Gram-Schmidt with pivoting: keep vectors whose residual norm
     * stays above a threshold after projecting out previously-accepted ones,
     * until m_needed orthonormal vectors are found. */
    int accepted = 0;
    for (int t = 0; t < cluster_size && accepted < m_needed; t++) {
      complex_t *v = &z[t * n];

      for (int a = 0; a < accepted; a++) {
        complex_t dot = c_zero();

        for (int i = 0; i < n; i++) {
          dot =
              c_add(dot, c_mul(c_conj(CMAT(result->eigenvectors, i, k_out + a)),
                               v[i]));
        }

        for (int i = 0; i < n; i++) {
          v[i] =
              c_sub(v[i], c_mul(dot, CMAT(result->eigenvectors, i, k_out + a)));
        }
      }

      double norm2 = 0.0;
      for (int i = 0; i < n; i++) {
        norm2 += c_abs2(v[i]);
      }

      double norm = sqrt(norm2);
      if (norm < 1e-10) {
        continue; // linearly dependent on already-accepted vectors
      }

      for (int i = 0; i < n; i++) {
        CMAT(result->eigenvectors, i, k_out + accepted) =
            c_scale(v[i], 1.0 / norm);
      }

      accepted++;
    }

    free(z);

    for (int a = 0; a < m_needed; a++) {
      result->eigenvalues[k_out + a] = lambda;
    }

    k_out += m_needed;
    pidx = group_end;
  }

  free(pair_lambda);
  free(pair_col0);
  free(pair_col1);
  free(order);
  eigen_free(eig2n);

  return result;

#endif
}
