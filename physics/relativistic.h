#ifndef QMC_RELATIVISTIC_H
#define QMC_RELATIVISTIC_H

#include "../core/matrix.h"
#include "../core/vector.h"

/* Solve 1D Klein-Gordon equation for scalar potential V(x):
   ( -\hbar^2 c^2 d^2/dx^2 + m^2 c^4 ) psi = (E - V)^2 \psi
   (after factoring out time dependence)
*/
eigen_t *klein_gordon_1d(double *x, int N, double *V, double m, double hbar,
                         double c);

/* Solve 1D Dirac equation (2-component) for scalar potential:
   [ c \sigma_z p + m c^2 \sigma_x + V(x) ] \phi = E \phi
   where \sigma_x, \sigma_z are Pauli matrices.
   Returns eigen_t with eigenvalues (energies) and eigenvectors (spinors).
*/
eigen_t *dirac_1d(double *x, int N, double *V, double m, double hbar, double c);

// TODO: Implement dirac_2d and dirac_3d

#endif
