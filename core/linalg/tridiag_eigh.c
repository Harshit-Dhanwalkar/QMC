/*
 * Tridiagonal symmetric eigensolver:
 * Implicit QL algorithm with Wilkinson shift (EISPACK tql2 / Numerical Recipes
 * "tqli", adapted to 0-indexed C).
 */

#include "tridiag_eigh.h"
#include "../complex.h"
#include "../matrix.h"
#include <math.h>
#include <stdlib.h>

static double pythag(double a, double b) {
  double absa = fabs(a), absb = fabs(b);
  if (absa > absb) {
    double r = absb / absa;
    return absa * sqrt(1.0 + r * r);
  } else if (absb == 0.0) {
    return 0.0;
  } else {
    double r = absa / absb;
    return absb * sqrt(1.0 + r * r);
  }
}

/*
 * Shared QL sweep. Computes ascending eigenvalues into returned array.
 * The per-rotation update loop below a stride-1 sweep over row, rather than
 * stride-n sweep a row-major layout would give). If z is NULL, all eigenvector
 * bookkeeping is skipped entirely.
 *
 * Returns NULL on allocation failure.
 */
static double *tridiag_eigh_core(const double *diag, const double *offdiag,
                                 int n, double *z) {
  double *d = malloc((size_t)n * sizeof *d);
  double *e = malloc((size_t)n * sizeof *e);
  if (!d || !e) {
    free(d);
    free(e);
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    d[i] = diag[i];
  }

  for (int i = 0; i < n - 1; i++) {
    e[i] = offdiag[i];
  }
  e[n - 1] = 0.0;

  for (int l = 0; l < n; l++) {
    int iter = 0;
    int m;
    do {
      for (m = l; m < n - 1; m++) {
        double dd = fabs(d[m]) + fabs(d[m + 1]);
        if (fabs(e[m]) <= dd * 1e-15)

          break;
      }

      if (m != l) {
        if (++iter > 100) {
          break; // did not converge for eigenvalue
        }

        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        double r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? fabs(r) : -fabs(r)));
        double s = 1.0, c = 1.0, p = 0.0;
        int i;
        for (i = m - 1; i >= l; i--) {
          double f = s * e[i];
          double b = c * e[i];
          r = pythag(f, g);
          e[i + 1] = r;
          if (r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;

            break;
          }

          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;

          if (z) {
            // Column-major (z[col*n+row]): sweep over row (k) is stride-1, not
            // stride-n
            for (int k = 0; k < n; k++) {
              double f2 = z[(i + 1) * n + k];
              z[(i + 1) * n + k] = s * z[i * n + k] + c * f2;
              z[i * n + k] = c * z[i * n + k] - s * f2;
            }
          }
        }

        if (r == 0.0 && i >= l) {
          continue;
        }

        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }

  // Sort ascending, permuting eigenvector columns to match (only if
  // eigenvectors were requested).
  for (int i = 0; i < n - 1; i++) {
    int k = i;
    double p = d[i];
    for (int j = i + 1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }

    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      if (z) {
        // Column-major: swap two contiguous n-length column blocks directly
        // (also a cache-friendly sequential sweep).
        for (int row = 0; row < n; row++) {
          double tmp = z[i * n + row];
          z[i * n + row] = z[k * n + row];
          z[k * n + row] = tmp;
        }
      }
    }
  }

  free(e);

  return d;
}

eigen_t *tridiag_eigh(const double *diag, const double *offdiag, int n) {
  if (!diag || n < 1 || (n > 1 && !offdiag)) {
    return NULL;
  }

  double *z = malloc((size_t)n * n * sizeof *z); // column-major
  if (!z) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      z[i * n + j] = (i == j) ? 1.0 : 0.0; // identity is symmetric either way
    }
  }

  double *d = tridiag_eigh_core(diag, offdiag, n, z);
  if (!d) {
    free(z);

    return NULL;
  }

  cmatrix_t *Z = cmatrix_alloc(n, n);
  if (!Z) {
    free(d);
    free(z);

    return NULL;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(Z, i, j) = c_real(z[j * n + i]); // column-major z -> row-major CMAT
    }
  }
  free(z);

  eigen_t *result = malloc(sizeof *result);
  if (!result) {
    free(d);
    cmatrix_free(Z);

    return NULL;
  }

  result->n = n;
  result->eigenvalues = d;
  result->eigenvectors = Z;

  return result;
}

eigen_t *tridiag_eigvals(const double *diag, const double *offdiag, int n) {
  if (!diag || n < 1 || (n > 1 && !offdiag)) {
    return NULL;
  }

  double *d = tridiag_eigh_core(diag, offdiag, n, NULL);
  if (!d) {
    return NULL;
  }

  eigen_t *result = malloc(sizeof *result);
  if (!result) {
    free(d);

    return NULL;
  }

  result->n = n;
  result->eigenvalues = d;
  result->eigenvectors = NULL; // no eigenvectors were computed

  return result;
}
