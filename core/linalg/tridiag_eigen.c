/*
 * QR algorithm for symmetric tridiagonal matrices
 *
 * For full dense matrices (3D, coupled systems):
 * TODO:
 * - Householder tridiagonalization :
 *   https://en.wikipedia.org/wiki/Householder_transformation
 * - QR iteration (Golub & Van Loan algorithm) :
 *   https://github.com/birocoles/matcomp
 *
 * HACK: Use LAPACK via eigen_generic.
 */

#include "../complex.h"
#include "../matrix.h"
#include "../vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* TODO:
// QR decomposition using Householder reflections
// HACK: Reduces Hermitian matrix to tridiagonal form first for efficiency
static void householder_reflection(cmatrix_t *A, int k) {
  int m = A->nrows - k;
  if (m <= 1) {
    return;
  }

  // Extract column k below diagonal
  complex_t *x = malloc(m * sizeof(complex_t));
  for (int i = 0; i < m; i++) {
    x[i] = CMAT(A, k + i, k);
  }

  // Compute Householder vector
  double norm = 0.0;
  for (int i = 0; i < m; i++) {
    norm += c_abs2(x[i]);
  }
  norm = sqrt(norm);

  if (norm < 1e-15) {
    free(x);

    return;
  }

  complex_t sigma = c_mul(x[0], c_real(norm / c_abs(x[0])));
  x[0] = c_add(x[0], sigma);

  double v_norm = 0.0;
  for (int i = 0; i < m; i++) {
    v_norm += c_abs2(x[i]);
  }
  v_norm = sqrt(v_norm);

  if (v_norm < 1e-15) {
    free(x);

    return;
  }

  // Apply reflection: A = (I - 2vv^\dagger / ||v||^2)A(I - 2vv^\dagger /
  // ||v||^2)
  for (int i = 0; i < m; i++) {
    x[i] = c_scale(x[i], 1.0 / v_norm);
  }

  // Left multiplication
  for (int j = k; j < A->ncols; j++) {
    complex_t sum = c_zero();

    for (int i = 0; i < m; i++) {
      sum = c_add(sum, c_mul(c_conj(x[i]), CMAT(A, k + i, j)));
    }
    complex_t factor = c_scale(sum, -2.0);

    for (int i = 0; i < m; i++) {
      CMAT(A, k + i, j) = c_add(CMAT(A, k + i, j), c_mul(x[i], factor));
    }
  }

  // Right multiplication
  for (int i = 0; i < A->nrows; i++) {
    complex_t sum = c_zero();

    for (int j = 0; j < m; j++) {
      sum = c_add(sum, c_mul(CMAT(A, i, k + j), c_conj(x[j])));
    }
    complex_t factor = c_scale(sum, -2.0);

    for (int j = 0; j < m; j++) {
      CMAT(A, i, k + j) = c_add(CMAT(A, i, k + j), c_mul(factor, x[j]));
    }
  }

  free(x);
}
*/

/*  TODO:
// Power iteration for lowest eigenvalue/eigenvector
static void power_iteration(cmatrix_t *A, cvector_t *v, double *lambda,
                            int max_iter, double tol) {
  cmatrix_t *Acopy = cmatrix_copy(A);
  cvector_t *w = cvector_alloc(A->nrows);
  cvector_t *v_old = cvector_copy(v);

  for (int iter = 0; iter < max_iter; iter++) {
    // w = A @ v
    for (int i = 0; i < A->nrows; i++) {
      w->data[i] = c_zero();

      for (int j = 0; j < A->ncols; j++) {
        w->data[i] = c_add(w->data[i], c_mul(CMAT(A, i, j), v->data[j]));
      }
    }

    // Rayleigh quotient: \lambda = v^\dagger A v / (v^\dagger v)
    complex_t numerator = c_zero();
    complex_t denominator = c_zero();
    for (int i = 0; i < A->nrows; i++) {
      numerator = c_add(numerator, c_mul(c_conj(v->data[i]), w->data[i]));
      denominator = c_add(denominator, c_mul(c_conj(v->data[i]), v->data[i]));
    }
    *lambda = c_div(numerator, denominator).re; // real part is the eigenvalue

    // Normalize w
    cvector_normalize(w);

    // Convergence check
    double err = 0.0;
    for (int i = 0; i < A->nrows; i++) {
      complex_t diff = c_sub(w->data[i], v_old->data[i]);
      err += c_abs2(diff);
    }

    err = sqrt(err);
    if (err < tol) {
      break;
    }

    memcpy(v_old->data, v->data, A->nrows * sizeof(complex_t));
    memcpy(v->data, w->data, A->nrows * sizeof(complex_t));
  }

  cmatrix_free(Acopy);
  cvector_free(w);
  cvector_free(v_old);
}
*/

/* TODO:
// Inverse iteration for specific eigenvalue
static void inverse_iteration(cmatrix_t *A, cvector_t *v, double sigma,
int max_iter, double tol) {
cmatrix_t *B = cmatrix_copy(A);

// Shift: B = A - \sigma I
for (int i = 0; i < B->nrows; i++) {
  CMAT(B, i, i) = c_sub(CMAT(B, i, i), c_real(sigma));
}

cvector_t *w = cvector_alloc(A->nrows);
cvector_t *v_old = cvector_copy(v);

for (int iter = 0; iter < max_iter; iter++) {
  // Solve (A - \sigma I)y = v
  cvector_t *y = cmatrix_solve(B, v);
  if (!y) {
    break;
  }

  cvector_normalize(y);

  // Convergence check
  double err = 0.0;
  for (int i = 0; i < A->nrows; i++) {
    complex_t diff = c_sub(y->data[i], v_old->data[i]);
    err += c_abs2(diff);
  }
  err = sqrt(err);

  if (err < tol) {
    memcpy(v->data, y->data, A->nrows * sizeof(complex_t));
    cvector_free(y);

    break;
  }

  memcpy(v_old->data, v->data, A->nrows * sizeof(complex_t));
  memcpy(v->data, y->data, A->nrows * sizeof(complex_t));

  cvector_free(y);
}

cmatrix_free(B);
cvector_free(w);
cvector_free(v_old);
}
*/

// Jacobi eigenvalue algorithm for real symmetric matrices
static double max_offdiag(const cmatrix_t *A) {
  double max = 0.0;
  int n = A->nrows;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double val = fabs(CMAT(A, i, j).re);

      if (val > max) {
        max = val;
      }
    }
  }

  return max;
}

// Sum of squares of all off-diagonal elements
static double offdiag_norm(const cmatrix_t *A) {
  double sum = 0.0;
  int n = A->nrows;

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      sum += CMAT(A, i, j).re * CMAT(A, i, j).re;
    }
  }

  return sqrt(2.0 * sum);
}

static void jacobi_rotate(cmatrix_t *A, cmatrix_t *V, int i, int j) {
  int n = A->nrows;
  double a_ii = CMAT(A, i, i).re;
  double a_jj = CMAT(A, j, j).re;
  double a_ij = CMAT(A, i, j).re;

  if (fabs(a_ij) < 1e-300) {
    return; // already zero
  }

  double theta = 0.5 * atan2(2.0 * a_ij, a_jj - a_ii);
  double c = cos(theta);
  double s = sin(theta);

  // Apply rotation to A (similarity transform)
  // Update rows/cols i and j for all k != i,j
  for (int k = 0; k < n; k++) {
    if (k != i && k != j) {
      double a_ki = CMAT(A, k, i).re;
      double a_kj = CMAT(A, k, j).re;
      double new_ki = c * a_ki - s * a_kj;
      double new_kj = s * a_ki + c * a_kj;
      CMAT(A, k, i) = c_real(new_ki);
      CMAT(A, i, k) = c_real(new_ki);
      CMAT(A, k, j) = c_real(new_kj);
      CMAT(A, j, k) = c_real(new_kj);
    }
  }

  // Update diagonal and off-diagonal of i,j
  double new_ii = c * c * a_ii - 2.0 * s * c * a_ij + s * s * a_jj;
  double new_jj = s * s * a_ii + 2.0 * s * c * a_ij + c * c * a_jj;
  CMAT(A, i, i) = c_real(new_ii);
  CMAT(A, j, j) = c_real(new_jj);
  CMAT(A, i, j) = c_real(0.0);
  CMAT(A, j, i) = c_real(0.0);

  // Accumulate rotation into eigenvector matrix: V = V * R
  for (int k = 0; k < n; k++) {
    double v_ki = CMAT(V, k, i).re;
    double v_kj = CMAT(V, k, j).re;
    CMAT(V, k, i) = c_real(c * v_ki - s * v_kj);
    CMAT(V, k, j) = c_real(s * v_ki + c * v_kj);
  }
}

// Public API
// Main eigendecomposition for Hermitian matrices
eigen_t *cmatrix_eigh(const cmatrix_t *A) {
  if (!A || A->nrows != A->ncols) {
    fprintf(stderr, "Error: matrix must be square\n");
    // fprintf(stderr, "Error: matrix must be square for eigen
    // decomposition\n");

    return NULL;
  }

  int n = A->nrows;

  /*
  eigen_t *result = malloc(sizeof(eigen_t));
  if (!result) {
    return NULL;
  }

  result->n = n;
  result->eigenvalues = malloc(n * sizeof(double));
  result->eigenvectors = cmatrix_alloc(n, n); // allocate an nxn matrix
  if (!result->eigenvalues || !result->eigenvectors) {
    free(result->eigenvalues);
    free(result->eigenvectors);
    free(result);

    return NULL;
  }

  // Compute all eigenvalues/eigenvectors
  // TODO: use LAPACK/BLAS
  cmatrix_t *Acopy = cmatrix_copy(A);
  if (!Acopy) {
    eigen_free(result);
    cmatrix_free(Acopy);

    return NULL;
  }
  */

  /* FIX:
  // Power iteration for ground state
  cvector_t *v = cvector_alloc(n);
  for (int i = 0; i < n; i++) {
    v->data[i] = c_real(1.0 / sqrt(n)); // Random initial guess
  }

  power_iteration(Acopy, v, &result->eigenvalues[0], 1000, 1e-10);
  result->eigenvectors[0] = *v;
  free(v);

  // TODO: For higher states: use Ritz values or Lanczos method
  // HACK: For now, approximate with shifted power iteration
  for (int k = 1; k < (n < 5 ? n : 5); k++) {
    v = cvector_alloc(n);

    // Random orthogonal starting vector
    for (int i = 0; i < n; i++) {
      v->data[i] = c_real((double)rand() / RAND_MAX);
    }

    // Gram-Schmidt orthogonalize against previous eigenvectors
    for (int j = 0; j < k; j++) {
      complex_t proj = c_zero();
      for (int i = 0; i < n; i++) {
        proj = c_add(
            proj, c_mul(c_conj(result->eigenvectors[j].data[i]), v->data[i]));
      }

      for (int i = 0; i < n; i++) {
        v->data[i] =
            c_sub(v->data[i], c_mul(proj, result->eigenvectors[j].data[i]));
      }
    }
    cvector_normalize(v);

    power_iteration(Acopy, v, &result->eigenvalues[k], 1000, 1e-10);
    result->eigenvectors[k] = *v;

    free(v);
  }
  */

  /* HACK: Compute up to min(n, 5) eigenpairs
  int num_states = (n < 5) ? n : 5;

  for (int k = 0; k < num_states; k++) {
    cvector_t *v = cvector_alloc(n);
    if (!v) {
      eigen_free(result);
      cmatrix_free(Acopy);

      return NULL;
    }

    if (k == 0) {
      // Ground state: uniform initial guess
      for (int i = 0; i < n; i++) {
        v->data[i] = c_real(1.0 / sqrt(n));
      }
    } else {
      // Higher states: random vector, orthogonalized against previous
      for (int i = 0; i < n; i++) {
        v->data[i] = c_real((double)rand() / RAND_MAX);
      }

      for (int j = 0; j < k; j++) {
        complex_t proj = c_zero();

        for (int i = 0; i < n; i++) {
          proj = c_add(proj, c_mul(c_conj(CMAT(result->eigenvectors, i, j)),
                                   v->data[i]));
        }

        for (int i = 0; i < n; i++) {
          v->data[i] =
              c_sub(v->data[i], c_mul(proj, CMAT(result->eigenvectors, i, j)));
        }
      }

      cvector_normalize(v);
    }

    double lambda;
    power_iteration(Acopy, v, &lambda, 1000, 1e-10);
    result->eigenvalues[k] = lambda;

    // Store eigenvector as column k of the matrix
    for (int i = 0; i < n; i++) {
      CMAT(result->eigenvectors, i, k) = v->data[i];
    }

    cvector_free(v);
  }

  // For remaining eigenvalues, set to zero
  for (int k = num_states; k < n; k++) {
    // result->eigenvalues[k] = 0.0;
    for (int i = 0; i < n; i++) {
      CMAT(result->eigenvectors, i, k) = c_zero();
    }
  }
  */

  // Copy input to working matrix
  cmatrix_t *Acopy = cmatrix_copy(A);
  if (!Acopy) {
    return NULL;
  }

  // Initialize eigenvector matrix to identity
  cmatrix_t *V = cmatrix_alloc(n, n);
  if (!V) {
    cmatrix_free(Acopy);

    return NULL;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(V, i, j) = (i == j) ? c_real(1.0) : c_real(0.0);
    }
  }

  /* Run the Jacobi sweep until convergence.
   * Cyclic Jacobi: sweep all (p,q) pairs with p<q, repeat until
   * off-diagonal norm is below tolerance.
   */
  const double tol = 1e-12;
  const int max_sweeps = 100;
  int sweep = 0;

  while (offdiag_norm(Acopy) > tol * (double)n && sweep < max_sweeps) {
    for (int p = 0; p < n - 1; p++) {
      for (int q = p + 1; q < n; q++) {
        jacobi_rotate(Acopy, V, p, q);
      }
    }

    sweep++;
  }

  if (sweep >= max_sweeps) {
    fprintf(stderr,
            "Warning cmatrix_eigh: did not fully converge after %d sweeps "
            "off-diag norm = %.3e).\n Result may be approximate.\n",
            max_sweeps, offdiag_norm(Acopy));
  }

  // Build result
  eigen_t *result = malloc(sizeof *result);
  if (!result) {
    cmatrix_free(Acopy);
    cmatrix_free(V);

    return NULL;
  }

  result->n = n;
  result->eigenvalues = malloc(n * sizeof(double));
  result->eigenvectors = V;
  if (!result->eigenvalues) {
    free(result);
    cmatrix_free(Acopy);
    cmatrix_free(V);

    return NULL;
  }

  // Copy diagonal as eigenvalues (sorted in descending order by Jacobi)
  for (int i = 0; i < n; i++) {
    result->eigenvalues[i] = CMAT(Acopy, i, i).re;
  }

  // Sort eigenvalues and eigenvectors in ascending order
  // NOTE: n is of order 100
  for (int i = 0; i < n - 1; i++) {
    int min_idx = i;

    for (int j = i + 1; j < n; j++) {
      if (result->eigenvalues[j] < result->eigenvalues[min_idx]) {
        min_idx = j;
      }
    }

    if (min_idx != i) {
      double tmp = result->eigenvalues[i];

      result->eigenvalues[i] = result->eigenvalues[min_idx];
      result->eigenvalues[min_idx] = tmp;
      for (int k = 0; k < n; k++) {
        complex_t ctmp = CMAT(V, k, i);
        CMAT(V, k, i) = CMAT(V, k, min_idx);
        CMAT(V, k, min_idx) = ctmp;
      }
    }
  }

  cmatrix_free(Acopy);

  return result;
}

void eigen_free(eigen_t *e) {
  if (!e) {
    return;
  }

  free(e->eigenvalues);
  cmatrix_free(e->eigenvectors);
  free(e);
}
