#ifndef QMC_TRIDIAG_EIGH_H
#define QMC_TRIDIAG_EIGH_H

#include "../matrix.h"

/*
 * Eigen-decomposition of a real symmetric TRIDIAGONAL matrix via the
 * implicit QL algorithm with Wilkinson shift (EISPACK tql2 /
 * Numerical-Recipes "tqli").
 *
 * diag:     diagonal entries, size n
 * offdiag:  off-diagonal entries, size n-1;
 *           offdiag[i] connects row i and row i+1
 *
 * Returns eigen_t* with eigenvalues ascending and eigenvectors as
 * columns of an nxn cmatrix_t.
 */
eigen_t *tridiag_eigh(const double *diag, const double *offdiag, int n);

/*
 * Eigenvalues only - same as tridiag_eigh, but skips all eigenvector
 * bookkeeping entirely (no working-array allocation, no per-rotation update, no
 * eigenvector permutation during sorting).
 *
 * Returns eigen_t* with eigenvalues ascending and eigenvectors == NULL
 * Used instead of tridiag_eigh whenever caller only reads eig->eigenvalues.
 */
eigen_t *tridiag_eigvals(const double *diag, const double *offdiag, int n);

#endif
