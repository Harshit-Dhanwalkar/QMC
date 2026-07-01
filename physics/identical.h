#ifndef QMC_IDENTICAL_H
#define QMC_IDENTICAL_H

#include "../core/matrix.h"
#include "../core/vector.h"

/* Slater determinant for fermions: given N single-particle states (as vectors),
   returns a matrix where each column is orbital and determinant is N-particle
   wavefunction. This is just matrix of orbitals.
*/
cmatrix_t *slater_determinant(cvector_t **orbitals, int N);

/* Symmetrize product of orbitals for bosons (permanent) */
cmatrix_t *bosonic_symmetrize(cvector_t **orbitals, int N);

#endif
