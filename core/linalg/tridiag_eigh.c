/*
 * Tridiagonal symmetric eigensolver: implicit QL algorithm with
 * Wilkinson shift (EISPACK tql2 / Numerical Recipes "tqli",

For full dense matrices (3D, coupled systems):
- TODO: Householder tridiagonalization -- Commented out
(https://en.wikipedia.org/wiki/Householder_transformation)
- QR iteration (Golub & Van Loan algorithm :
https://github.com/birocoles/matcomp)

// HACK: Use LAPACK via eigen_generic.
*/

#include "../complex.h"
#include "../matrix.h"
#include "tridiag_eigh.h"
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

eigen_t *tridiag_eigh(const double *diag, const double *offdiag, int n) {
  if (!diag || n < 1 || (n > 1 && !offdiag))
    return NULL;

  double *d = malloc((size_t)n * sizeof *d); // eigenvalues on exit
  double *e = malloc((size_t)n * sizeof *e); // QL working array
  if (!d || !e) {
    free(d);
    free(e);
    return NULL;
  }
  for (int i = 0; i < n; i++)
    d[i] = diag[i];
  // Build e[] in the QL sweep's convention
  e[0] = 0.0;
  for (int i = 1; i < n; i++)
    e[i] = offdiag[i - 1];
  for (int i = 1; i < n; i++)
    e[i - 1] = e[i];
  e[n - 1] = 0.0;

  double *z = malloc((size_t)n * n * sizeof *z); // row-major eigenvectors
  if (!z) {
    free(d);
    free(e);
    return NULL;
  }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      z[i * n + j] = (i == j) ? 1.0 : 0.0;

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
        if (++iter > 100)
          break; // did not converge for eigenvalue

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
          for (int k = 0; k < n; k++) {
            double f2 = z[k * n + (i + 1)];
            z[k * n + (i + 1)] = s * z[k * n + i] + c * f2;
            z[k * n + i] = c * z[k * n + i] - s * f2;
          }
        }
        if (r == 0.0 && i >= l)
          continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }

  // Sort ascending, permuting eigenvector columns to match
  for (int i = 0; i < n - 1; i++) {
    int k = i;
    double p = d[i];
    for (int j = i + 1; j < n; j++)
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < n; j++) {
        double tmp = z[j * n + i];
        z[j * n + i] = z[j * n + k];
        z[j * n + k] = tmp;
      }
    }
  }

  cmatrix_t *Z = cmatrix_alloc(n, n);
  if (!Z) {
    free(d);
    free(e);
    free(z);
    return NULL;
  }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      CMAT(Z, i, j) = c_real(z[i * n + j]);
  free(z);

  eigen_t *result = malloc(sizeof *result);
  if (!result) {
    free(d);
    free(e);
    cmatrix_free(Z);
    return NULL;
  }
  result->n = n;
  result->eigenvalues = d;
  result->eigenvectors = Z;

  free(e);
  return result;
}
