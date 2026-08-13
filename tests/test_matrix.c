#include "../core/complex.h"
#include "../core/matrix.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
  printf(" > Testing matrix operations...\n");

  // Allocate a 2x3 matrix
  cmatrix_t *A = cmatrix_alloc(2, 3);
  if (!A) {
    printf("FAIL: allocation\n");
    return 1;
  }

  CMAT(A, 0, 0) = c_real(1.0);
  CMAT(A, 0, 1) = c_real(2.0);
  CMAT(A, 0, 2) = c_real(3.0);
  CMAT(A, 1, 0) = c_real(4.0);
  CMAT(A, 1, 1) = c_real(5.0);
  CMAT(A, 1, 2) = c_real(6.0);

  // Test transpose (2x3 -> 3x2)
  cmatrix_t *At = cmatrix_transpose(A);
  if (!At) {
    printf("FAIL: transpose\n");

    cmatrix_free(A);

    return 1;
  }
  if (At->nrows != 3 || At->ncols != 2) {
    printf("FAIL: transpose dimensions\n");

    cmatrix_free(A);
    cmatrix_free(At);

    return 1;
  }
  if (c_abs(c_sub(CMAT(At, 0, 0), c_real(1.0))) > 1e-12) {
    printf("FAIL: transpose value\n");

    cmatrix_free(A);
    cmatrix_free(At);

    return 1;
  }

  // Test multiplication (A * At) -> 2x2
  cmatrix_t *B = cmatrix_multiply(A, At);
  if (!B) {
    printf("FAIL: multiplication\n");

    cmatrix_free(A);
    cmatrix_free(At);

    return 1;
  }

  if (B->nrows != 2 || B->ncols != 2) {
    printf("FAIL: multiplication dimensions\n");

    cmatrix_free(A);
    cmatrix_free(At);
    cmatrix_free(B);

    return 1;
  }

  // Expected: [[14, 32], [32, 77]]
  if (c_abs(c_sub(CMAT(B, 0, 0), c_real(14.0))) > 1e-12) {
    printf("FAIL: multiplication result\n");

    cmatrix_free(A);
    cmatrix_free(At);
    cmatrix_free(B);

    return 1;
  }

  // Free everything
  cmatrix_free(A);
  cmatrix_free(At);
  cmatrix_free(B);

  return 0;
}
