#ifndef QMC_ANGULAR_H
#define QMC_ANGULAR_H

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"

// Pauli matrices (2x2)
extern complex_t sigma_x[4]; // row-major, indexed as [i*2 + j]
extern complex_t sigma_y[4];
extern complex_t sigma_z[4];

// Ladder operators for angular momentum: L_{\pm} = L_x \pm i L_y
complex_t l_plus_op(int l, int m,
                    int m_prime); // matrix element <l,m+1|L_+|l,m>
complex_t l_minus_op(int l, int m, int m_prime); // <l,m-1|L_-|l,m>

// Matrix representation of L^2, L_z for a given l (dimension 2l+1)
cmatrix_t *l2_matrix(int l);
cmatrix_t *lz_matrix(int l);
cmatrix_t *lx_matrix(int l);
cmatrix_t *ly_matrix(int l);

// Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M> for small j (j1,j2 <=
//   Returns real coefficient (phase convention Condon-Shortley)
double clebsch_gordan(int j1, int m1, int j2, int m2, int J, int M);

// Spin-1/2 operations: apply sigma_x, sigma_y, sigma_z to a 2-component spi
void spin_sigma_x(cvector_t *spinor);
void spin_sigma_y(cvector_t *spinor);
void spin_sigma_z(cvector_t *spinor);

#endif
