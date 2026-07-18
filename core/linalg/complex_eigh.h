#ifndef QMC_COMPLEX_EIGH_H
#define QMC_COMPLEX_EIGH_H

#include "../matrix.h"

/*
 * Eigen-decomposition of general complex Hermitian matrix via real-embedding
 * trick.
 *
 * For H = A + iB, build real symmetric 2N x 2N matrix
 *   M = [[ A, -B ],
 *        [ B,  A ]]
 *
 */
eigen_t *cmatrix_eigh_complex(cmatrix_t *H);

#endif
