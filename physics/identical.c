/*
Identical particles: Slater determinants and symmetrization.
*/

#include "identical.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

cmatrix_t *slater_matrix(cvector_t **orbitals, int N, const int *indices) {
  if (!orbitals || !indices || N < 1) {
    return NULL;
  }

  cmatrix_t *M = cmatrix_alloc(N, N);
  if (!M) {
    return NULL;
  }

  for (int i = 0; i < N; i++) {
    if (!orbitals[i] || indices[i] < 0 || indices[i] >= orbitals[i]->n) {
      cmatrix_free(M);

      return NULL;
    }
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CMAT(M, i, j) = orbitals[i]->data[indices[j]];
    }
  }

  return M;
}

// Complex determinant via Gaussian elimination with partial pivoting
static complex_t determinant_gauss_destructive(cmatrix_t *M) {
  int n = M->nrows;
  complex_t det = c_real(1.0);

  for (int col = 0; col < n; col++) {
    int pivot = col;
    double best = c_abs2(CMAT(M, col, col));
    for (int r = col + 1; r < n; r++) {
      double mag = c_abs2(CMAT(M, r, col));
      if (mag > best) {
        best = mag;
        pivot = r;
      }
    }

    if (best < 1e-300) {
      return c_zero(); // singular matrix
    }

    if (pivot != col) {
      for (int c = 0; c < n; c++) {
        complex_t tmp = CMAT(M, col, c);
        CMAT(M, col, c) = CMAT(M, pivot, c);
        CMAT(M, pivot, c) = tmp;
      }

      det = c_scale(det, -1.0); // row swap flips determinant's sign
    }

    complex_t piv = CMAT(M, col, col);
    det = c_mul(det, piv);
    complex_t piv_recip = c_scale(c_conj(piv), 1.0 / c_abs2(piv));

    for (int r = col + 1; r < n; r++) {
      complex_t factor = c_mul(CMAT(M, r, col), piv_recip);
      for (int c = col; c < n; c++) {
        CMAT(M, r, c) = c_sub(CMAT(M, r, c), c_mul(factor, CMAT(M, col, c)));
      }
    }
  }

  return det;
}

// Ryser's formula
// Complex permanent via Ryser's formula (Gray-code subset enumeration):
//
// NOTE: Complexity-theoretic basis of boson-sampling argument for quantum
// computational advantage (Aaronson & Arkhipov 2011)
static complex_t permanent_ryser(cmatrix_t *A) {
  int n = A->nrows;
  if (n == 0) {
    return c_real(1.0);
  }

  double *row_sum_re = calloc(n, sizeof *row_sum_re);
  double *row_sum_im = calloc(n, sizeof *row_sum_im);
  if (!row_sum_re || !row_sum_im) {
    free(row_sum_re);
    free(row_sum_im);

    return c_zero();
  }

  complex_t perm = c_zero();
  unsigned long long num_subsets = 1ULL << n;
  unsigned long long prev_gray = 0;

  for (unsigned long long k = 1; k < num_subsets; k++) {
    unsigned long long gray = k ^ (k >> 1);
    unsigned long long diff = gray ^ prev_gray;
    int col = 0;

    while (!((diff >> col) & 1ULL)) {
      col++;
    }

    int bit_turned_on = (gray >> col) & 1ULL;
    double s = bit_turned_on ? 1.0 : -1.0;

    for (int i = 0; i < n; i++) {
      row_sum_re[i] += s * CMAT(A, i, col).re;
      row_sum_im[i] += s * CMAT(A, i, col).im;
    }

    int subset_size = __builtin_popcountll(gray);
    double term_sign = (subset_size % 2 == 0) ? 1.0 : -1.0;

    complex_t prod = c_real(1.0);
    for (int i = 0; i < n; i++) {
      complex_t rs = {row_sum_re[i], row_sum_im[i]};
      prod = c_mul(prod, rs);
    }

    perm = c_add(perm, c_scale(prod, term_sign));
    prev_gray = gray;
  }

  free(row_sum_re);
  free(row_sum_im);

  double overall_sign = (n % 2 == 0) ? 1.0 : -1.0;

  return c_scale(perm, overall_sign);
}

static double factorial(int n) {
  double r = 1.0;
  for (int i = 2; i <= n; i++) {
    r *= i;
  }

  return r;
}

complex_t slater_determinant_value(cvector_t **orbitals, int N,
                                   const int *indices) {
  cmatrix_t *M = slater_matrix(orbitals, N, indices);
  if (!M) {
    return c_zero();
  }

  complex_t det = determinant_gauss_destructive(M);
  cmatrix_free(M);
  double inv_sqrt_nfact = 1.0 / sqrt(factorial(N));

  return c_scale(det, inv_sqrt_nfact);
}

complex_t bosonic_permanent_value(cvector_t **orbitals, int N,
                                  const int *indices) {
  cmatrix_t *M = slater_matrix(orbitals, N, indices);
  if (!M) {
    return c_zero();
  }

  complex_t perm = permanent_ryser(M);
  cmatrix_free(M);

  double mult_factor = 1.0;
  for (int i = 0; i < N; i++) {
    int already_counted = 0;
    for (int k = 0; k < i; k++) {
      if (orbitals[k] == orbitals[i]) {
        already_counted = 1;

        break;
      }
    }

    if (already_counted) {
      continue;
    }

    int count = 1;
    for (int k = i + 1; k < N; k++) {
      if (orbitals[k] == orbitals[i]) {
        count++;
      }
    }

    mult_factor *= factorial(count);
  }

  double inv_norm = 1.0 / sqrt(factorial(N) * mult_factor);

  return c_scale(perm, inv_norm);
}
