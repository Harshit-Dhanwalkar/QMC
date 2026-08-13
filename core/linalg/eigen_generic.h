#ifndef QMC_EIGEN_GENERIC_H
#define QMC_EIGEN_GENERIC_H

#include "../matrix.h"

/* Generic Hermitian eigen-solver.
   Uses LAPACK if available (USE_LAPACK), otherwise falls back to
   cmatrix_eigh (QR).
*/
eigen_t *cmatrix_eigh_generic(const cmatrix_t *A);

#endif
