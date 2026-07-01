#ifndef QMC_SVD_H
#define QMC_SVD_H

#include "../matrix.h"
#include "../vector.h"

/* SVD: A = U * diag(S) * V^\dagger.
   For complex A, we convert to real block matrix.
   U and V are complex matrices, S is real vector.
   Returns 0 on success.
*/
int svd_decompose(const cmatrix_t *A, cmatrix_t *U, cvector_t *S, cmatrix_t *V);

/* Solve least squares using SVD: min ||A x - b||.
   SVD must have been computed.
   Returns 0 on success.
*/
int svd_solve(const cmatrix_t *U, const cvector_t *S, const cmatrix_t *V,
              const cvector_t *b, cvector_t *x, double tol);

/* Pseudoinverse: A+ = V * diag(1/S) * U^\dagger.
   Returns new matrix.
*/
cmatrix_t *svd_pseudoinverse(const cmatrix_t *U, const cvector_t *S,
                             const cmatrix_t *V, double tol);

#endif
