#ifndef QMC_LU_H
#define QMC_LU_H

#include "../matrix.h"
#include "../vector.h"

/* LU decomposition with partial pivoting: PA = LU.
   Returns pivot array (int*) of size n, or NULL on failure.
   The input matrix A is overwritten with L+U (unit diagonal not stored).
   The pivot array indicates row permutations.
*/
int *lu_decompose(cmatrix_t *A);

/* Solve A x = b using LU decomposition.
   LU and pivot must come from lu_decompose.
   b is the right-hand side; x is the solution (must be preallocated).
   Returns 0 on success, -1 on failure.
*/
int lu_solve(const cmatrix_t *LU, const int *pivot, const cvector_t *b,
             cvector_t *x);

/* Determinant of A given LU and pivot.
   Returns complex determinant.
*/
complex_t lu_det(const cmatrix_t *LU, const int *pivot);

/* Inverse of A using LU decomposition.
   Returns new cmatrix_t* or NULL on failure.
*/
cmatrix_t *lu_invert(const cmatrix_t *A);

#endif
