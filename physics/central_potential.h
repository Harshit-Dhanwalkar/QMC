#ifndef QMC_CENTRAL_POTENTIAL_H
#define QMC_CENTRAL_POTENTIAL_H

#include "../core/linalg/tridiag_eigh.h"
#include "potentials.h"

/*
 * General 3D central-potential radial solver.
 *
 * Solves the radial Schrodinger equation for u(r) = r * R(r):
 *   -\hbar^2/(2m) u''(r) + [ V(r) + \hbar^2/(2m) * l(l+1)/r^2 ] u(r) = E u(r)
 *
 * on a uniform grid r[0..N-1], via finite-difference tridiagonal
 * discretization + tridiag_eigh (implicit QL).
 */
eigen_t *central_potential_radial_solve(double *r, int N, int l, double hbar,
                                        double mass, potential_fn V,
                                        void *params);

#endif // QMC_CENTRAL_POTENTIAL_H
