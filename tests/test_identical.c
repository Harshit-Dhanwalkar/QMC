/*
 * Gaussian-elimination determinant, Ryser's formula permanent
 *
 * 1. Random complex NxN matrices, N=2..6: determinant/permanent must
 *    match reference to numerical precision.
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/identical.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Not using Ryser formula for permanent
static complex_t cofactor_expand_reference(cmatrix_t *A, int use_sign) {
  int n = A->nrows;
  if (n == 1) {
    return CMAT(A, 0, 0);
  }

  complex_t sum = c_zero();
  double sign = 1.0;
  for (int j = 0; j < n; j++) {
    cmatrix_t *minor = cmatrix_alloc(n - 1, n - 1);

    for (int r = 1; r < n; r++) {
      int mc = 0;

      for (int c = 0; c < n; c++) {
        if (c == j) {
          continue;
        }

        CMAT(minor, r - 1, mc) = CMAT(A, r, c);
        mc++;
      }
    }

    complex_t sub = cofactor_expand_reference(minor, use_sign);
    complex_t term = c_mul(CMAT(A, 0, j), sub);
    if (use_sign) {
      term = c_scale(term, sign);
    }

    sum = c_add(sum, term);
    cmatrix_free(minor);

    sign = -sign;
  }

  return sum;
}

static int check_close_complex(complex_t got, complex_t expected, double tol,
                               const char *label) {
  double err =
      sqrt(pow(got.re - expected.re, 2) + pow(got.im - expected.im, 2));
  printf("  %s: got=(%.6f,%.6f) expected=(%.6f,%.6f) err=%.2e\n", label, got.re,
         got.im, expected.re, expected.im, err);

  return err > tol;
}

// Build orbitals such that slater_matrix(orbitals, N, indices) reproduces
// an arbitrary NxN complex matrix M exactly, with indices = {0,1,...,N-1}
static cvector_t **matrix_to_orbitals(cmatrix_t *M, int N) {
  cvector_t **orbitals = malloc(N * sizeof *orbitals);
  for (int i = 0; i < N; i++) {
    orbitals[i] = cvector_alloc(N);

    for (int j = 0; j < N; j++) {
      orbitals[i]->data[j] = CMAT(M, i, j);
    }
  }

  return orbitals;
}

static void free_orbitals(cvector_t **orbitals, int N) {
  for (int i = 0; i < N; i++) {
    cvector_free(orbitals[i]);
  }

  free(orbitals);
}

static double factorial_local(int n) {
  double r = 1.0;
  for (int i = 2; i <= n; i++) {
    r *= i;
  }

  return r;
}

static int test_against_reference(int N, unsigned int seed) {
  srand(seed);
  cmatrix_t *M = cmatrix_alloc(N, N);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double re = (double)rand() / RAND_MAX - 0.5;
      double im = (double)rand() / RAND_MAX - 0.5;
      CMAT(M, i, j) = c_add(c_real(re), c_imag(im));
    }
  }

  complex_t det_ref_raw = cofactor_expand_reference(M, 1);
  complex_t perm_ref_raw = cofactor_expand_reference(M, 0);
  double inv_sqrt_nfact = 1.0 / sqrt(factorial_local(N));
  complex_t det_ref = c_scale(det_ref_raw, inv_sqrt_nfact);
  complex_t perm_ref = c_scale(perm_ref_raw, inv_sqrt_nfact);

  cvector_t **orbitals = matrix_to_orbitals(M, N);
  cmatrix_free(M);
  int *indices = malloc(N * sizeof *indices);
  for (int i = 0; i < N; i++) {
    indices[i] = i;
  }

  complex_t det_new = slater_determinant_value(orbitals, N, indices);
  complex_t perm_new = bosonic_permanent_value(orbitals, N, indices);

  free_orbitals(orbitals, N);
  free(indices);

  char label1[32], label2[32];
  snprintf(label1, sizeof label1, "N=%d determinant", N);
  snprintf(label2, sizeof label2, "N=%d permanent", N);

  double tol = 1e-8;
  int fail = check_close_complex(det_new, det_ref, tol, label1);
  fail |= check_close_complex(perm_new, perm_ref, tol, label2);

  return fail;
}

static int test_diagonal_case(void) {
  int N = 4;
  cmatrix_t *M = cmatrix_alloc(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CMAT(M, i, j) = (i == j) ? c_real(2.0 + i) : c_zero();
    }
  }

  cvector_t **orbitals = matrix_to_orbitals(M, N);
  cmatrix_free(M);
  const int indices[4] = {0, 1, 2, 3};

  complex_t det = slater_determinant_value(orbitals, N, indices);
  complex_t perm = bosonic_permanent_value(orbitals, N, indices);
  free_orbitals(orbitals, N);

  // Diagonal entries are 2,3,4,5 -> product = 120. Both det and perm of
  // diagonal matrix equal product of diagonal.
  double product = 2.0 * 3.0 * 4.0 * 5.0;
  double inv_sqrt_nfact = 1.0 / sqrt(factorial_local(N));
  complex_t expected = c_real(product * inv_sqrt_nfact);

  int fail = check_close_complex(det, expected, 1e-8, "diagonal determinant");
  fail |= check_close_complex(perm, expected, 1e-8, "diagonal permanent");

  return fail;
}

int main(void) {
  int failed = 0;

  printf("Random matrices vs. O(N!) reference implementation:\n");
  for (int N = 2; N <= 6; N++) {
    failed += test_against_reference(N, 100 + N);
  }

  printf("Diagonal matrix (hand-checkable):\n");
  failed += test_diagonal_case();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
