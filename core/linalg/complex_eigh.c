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
  const double degeneracy_tol_floor = 1e-15;
  int k_out = 0;
  int idx = 0;

  while (idx < m2 && k_out < n) {
    int cluster_start = idx;
    int cluster_end = idx + 1;

    const int max_cluster_size = 8;
    while (cluster_end < m2 && cluster_end - cluster_start < max_cluster_size) {
      double lo = eig2n->eigenvalues[order[cluster_end - 1]];
      double hi = eig2n->eigenvalues[order[cluster_end]];
      double local_scale = fabs(lo) > fabs(hi) ? fabs(lo) : fabs(hi);
      double gap_tol = degeneracy_tol_rel * local_scale + degeneracy_tol_floor;
      if (fabs(hi - lo) >= gap_tol) {
        break;
      }

      cluster_end++;
    }

    if (cluster_end - cluster_start >= max_cluster_size) {
      /* Round down to an even boundary so the fallback below still pairs
       * cleanly; leave any odd remainder for the next cluster to pick up. */
      cluster_end = cluster_start + ((cluster_end - cluster_start) / 2) * 2;
    }

    int cluster_size = cluster_end - cluster_start; // always even (>=2)
    int m_needed = cluster_size / 2;

    double lambda_sum = 0.0;
    for (int t = cluster_start; t < cluster_end; t++) {
      lambda_sum += eig2n->eigenvalues[order[t]];
    }

    double lambda = lambda_sum / cluster_size;

    // Complex candidate vectors z_t = x_t + i*y_t for every real eigenvector in
    // this cluster
    complex_t *z = malloc((size_t)cluster_size * n * sizeof(complex_t));
    if (!z) {
      free(order);
      cmatrix_free(result->eigenvectors);
      free(result->eigenvalues);
      free(result);
      eigen_free(eig2n);

      return NULL;
    }

    for (int t = 0; t < cluster_size; t++) {
      int col = order[cluster_start + t];

      for (int i = 0; i < n; i++) {
        double x = CMAT(eig2n->eigenvectors, i, col).re;
        double y = CMAT(eig2n->eigenvectors, i + n, col).re;

        z[t * n + i] = c_add(c_real(x), c_imag(y));
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
    idx = cluster_end;
  }

  free(order);
  eigen_free(eig2n);

  return result;

#endif
}
