#ifndef QMC_IDENTICAL_H
#define QMC_IDENTICAL_H

#include "../core/matrix.h"
#include "../core/vector.h"

/*
 * Identical particles: build the N x N matrix needed to evaluate an
 * N-particle Slater determinant (fermions) or permanent (bosons) at a
 * specific configuration, and compute that determinant/permanent.
 *
 * orbitals[i]: tabulated single-particle orbital i, sampled on a common
 *              position grid of length M.
 * indices[j] : grid index of particle j's position (0 <= indices[j] < M).
 * Two particles at same grid index (indices[j1]==indices[j2]) or two identical
 * orbitals both make Slater determinant vanish exactly (Pauli exclusion)
 */

/* M_ij = orbitals[i]->data[indices[j]], size N x N */
cmatrix_t *slater_matrix(cvector_t **orbitals, int N, const int *indices);

/* Fermions: \psi(x_1,...,x_N) = (1 / \sqrt(N!)) * det[\phi_i(x_j)]
 * Exactly 0 if two particles share a position or two orbitals coincide.
 */
complex_t slater_determinant_value(cvector_t **orbitals, int N,
                                   const int *indices);

/* Bosons: \psi(x_1,...,x_N) = (1 / \sqrt(N! * prod_k n_k!)) * perm[\phi_i(x_j)]
 * is invariant under exchanging any two particles. n_k is the occupation number
 * of each distinct orbital among orbitals[0..n-1]; this reduces to simpler
 * 1 / sqrt(N!) exactly when every orbital is distinct, but prod_k n_k! factor
 * matters whenever multiple bosons share an orbital
 */
complex_t bosonic_permanent_value(cvector_t **orbitals, int N,
                                  const int *indices);

#endif
