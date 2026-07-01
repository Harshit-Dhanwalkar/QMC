#ifndef QMC_QR_H
#define QMC_QR_H

#include "../matrix.h"
#include "../vector.h"

/* QR decomposition: A = Q R, where Q is unitary, R upper triangular.
   A is overwritten with R (upper triangular part) and the Householder vectors
   (lower part). Q is output as a separate matrix (must be preallocated or
   NULL). Returns 0 on success, -1 on failure.
*/
int qr_decompose(cmatrix_t *A, cmatrix_t *Q);

/* Solve least squares: min ||A x - b|| using QR.
   A must have been QR decomposed (with qr_decompose).
   b is RHS, x is solution (preallocated).
   Returns 0 on success.
*/
int qr_solve(const cmatrix_t *A, const cmatrix_t *Q, const cvector_t *b,
             cvector_t *x);

/* Invert square matrix via QR (if non-singular).
   Returns new matrix or NULL.
*/
cmatrix_t *qr_invert(const cmatrix_t *A);

#endif
