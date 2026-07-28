#ifndef QMC_ANGULAR_H
#define QMC_ANGULAR_H

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"

/* Pauli matrices (2x2) */
extern complex_t sigma_x[4]; // row-major, indexed as [i*2 + j]
extern complex_t sigma_y[4];
extern complex_t sigma_z[4];

/* Ladder operators for angular momentum: L_{\pm} = L_x \pm i L_y */
complex_t l_plus_op(int l, int m,
                    int m_prime); // matrix element <l,m+1|L_+|l,m>
complex_t l_minus_op(int l, int m, int m_prime); // <l,m-1|L_-|l,m>

/* Matrix representation of L^2, L_z for a given l (dimension 2l+1) */
cmatrix_t *l2_matrix(int l);
cmatrix_t *lz_matrix(int l);
cmatrix_t *lx_matrix(int l);
cmatrix_t *ly_matrix(int l);

/*
 * Clebsch-Gordan coefficient <j1 m1; j2 m2 | J M>, general j1/j2/J via
 * Racah's formula
 *
 *   spin-1/2, m=+1/2    -> j=1, m=+1
 *   orbital   l=1, m=0  -> j=2, m=0
 *   total     J=3/2     -> J=3
 *
 * Returns 0 if inputs violate triangle inequality, m1+m2 != M, |m| > j for any
 * of three, or required integer/half-integer parity (j1+j2+J must be an
 * integer).
 * Phase convention: Condon-Shortley.
 */
double clebsch_gordan(int j1_2, int m1_2, int j2_2, int m2_2, int J_2, int M_2);

/*
 * Angular momentum coupling: |j1,j2; J,M> in uncoupled product basis {|j1,m1>
 * (x) |j2,m2>}
 */
cvector_t *couple_states(int j1_2, int j2_2, int J_2, int M_2);

/*
 * Enumerate allowed total-J values (doubled) for given j1_2, j2_2, i.e.
 * J_2 = |j1_2-j2_2|, |j1_2-j2_2|+2, ..., j1_2+j2_2.
 */
int couple_allowed_J(int j1_2, int j2_2, int *J2_out);

/*
 * Gaunt coefficient c^k(l,m; l',m') = \sqrt(4 * \pi / (2k+1)) * \int Y_lm^*
 * Y_kq Y_l'm' d\Omega, q = m - m'. Built from clebsch_gordan() above via
 *   c^k(l,m;l',m') = (-1)^(m-m') * sqrt((2l+1)/(2l'+1))
 *                    * <l0;k0|l'0> * <l,-m;k,q|l',-m'>
 * l, m, k, l', m' are integers, since orbital angular momenta are always
 * integers.
 *
 * Returns 0 if triangle inequality |l-l'| <= k <= l+l' is violated, if l+k+l'
 * is odd or if |m|>l or |m'|>l'.
 */
double gaunt_coefficient(int l, int m, int k, int lp, int mp);

// Spin-1/2 operations: apply \sigma_x, \sigma_y, \sigma_z to 2-component spinor
void spin_sigma_x(cvector_t *spinor);
void spin_sigma_y(cvector_t *spinor);
void spin_sigma_z(cvector_t *spinor);

#endif
