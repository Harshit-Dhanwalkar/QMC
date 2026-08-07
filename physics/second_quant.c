/*
Jordan-Wigner fermion-to-qubit mapping + an independent direct Fock-space
Hamiltonian builder for cross-validation
*/

#include "second_quant.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include <stdlib.h>

static const complex_t I2[4] = {{1, 0}, {0, 0}, {0, 0}, {1, 0}};
static const complex_t Z2[4] = {{1, 0}, {0, 0}, {0, 0}, {-1, 0}};
static const complex_t SIGMA_PLUS[4] = {
    {0, 0}, {0, 0}, {1, 0}, {0, 0}}; // |1><0|: row-major [00,01,10,11]
static const complex_t SIGMA_MINUS[4] = {{0, 0}, {1, 0}, {0, 0}, {0, 0}};

/* Kronecker product of two dense matrices:
 *  (A tensor B)[i*nb+k][j*nb+l] = A[i][j]*B[k][l]
 *
 * Returns a freshly allocated cmatrix_t
 */
static cmatrix_t *kron2(const cmatrix_t *A, const cmatrix_t *B) {
  int ra = A->nrows, ca = A->ncols, rb = B->nrows, cb = B->ncols;
  cmatrix_t *out = cmatrix_alloc(ra * rb, ca * cb);

  for (int i = 0; i < ra; i++) {
    for (int j = 0; j < ca; j++) {
      complex_t aij = CMAT(A, i, j);

      for (int k = 0; k < rb; k++) {
        for (int l = 0; l < cb; l++) {
          CMAT(out, i * rb + k, j * cb + l) = c_mul(aij, CMAT(B, k, l));
        }
      }
    }
  }

  return out;
}

// N-fold Kronecker product of a sequence of 2x2 gates (row-major [4] arrays),
// gates[0] (x) gates[1] (x) ... (x) gates[n-1]
static cmatrix_t *kron_chain_2x2(const complex_t *const *gates, int n) {
  cmatrix_t *result = cmatrix_alloc(2, 2);
  for (int k = 0; k < 4; k++) {
    result->data[k] = gates[0][k];
  }

  for (int idx = 1; idx < n; idx++) {
    cmatrix_t *next = cmatrix_alloc(2, 2);

    for (int k = 0; k < 4; k++) {
      next->data[k] = gates[idx][k];
    }

    cmatrix_t *combined = kron2(result, next);
    cmatrix_free(result);
    cmatrix_free(next);
    result = combined;
  }

  return result;
}

static cmatrix_t *jw_operator(int mode, int n_modes, const complex_t *local) {
  if (mode < 0 || mode >= n_modes || n_modes < 1) {
    return NULL;
  }

  const complex_t **gates = malloc((size_t)n_modes * sizeof *gates);
  for (int k = 0; k < mode; k++) {
    gates[k] = Z2;
  }

  gates[mode] = local;
  for (int k = mode + 1; k < n_modes; k++) {
    gates[k] = I2;
  }

  cmatrix_t *result = kron_chain_2x2(gates, n_modes);
  free(gates);

  return result;
}

cmatrix_t *jw_creation_operator(int mode, int n_modes) {
  return jw_operator(mode, n_modes, SIGMA_PLUS);
}

cmatrix_t *jw_annihilation_operator(int mode, int n_modes) {
  return jw_operator(mode, n_modes, SIGMA_MINUS);
}

// Number of set bits with index < mode (mode 0 = leftmost/MSB), used for
// fermionic anticommutation sign when creating/annihilating at `mode`
static int count_bits_before(int state, int mode, int n_modes) {
  int count = 0;
  for (int k = 0; k < mode; k++) {
    int bitpos = n_modes - 1 - k;
    if ((state >> bitpos) & 1) {
      count++;
    }
  }

  return count;
}

/* Direct (bit-manipulation) fermionic annihilation
 *
 * Returns resulting state index and writes fermionic sign to *sign_out, or
 * Returns -1 (Pauli exclusion / already empty) if mode isn't occupied.
 */
static int direct_annihilate(int state, int mode, int n_modes, int *sign_out) {
  int bitpos = n_modes - 1 - mode;
  if (!((state >> bitpos) & 1)) {
    return -1;
  }

  int n_before = count_bits_before(state, mode, n_modes);
  *sign_out = (n_before % 2 == 0) ? 1 : -1;

  return state & ~(1 << bitpos);
}

static int direct_create(int state, int mode, int n_modes, int *sign_out) {
  int bitpos = n_modes - 1 - mode;
  if ((state >> bitpos) & 1) {
    return -1; // already occupied
  }

  int n_before = count_bits_before(state, mode, n_modes);
  *sign_out = (n_before % 2 == 0) ? 1 : -1;

  return state | (1 << bitpos);
}

cmatrix_t *second_quant_build_hopping_hamiltonian(int n_modes,
                                                  const double *epsilon,
                                                  double t, double U) {
  if (n_modes < 1 || !epsilon) {
    return NULL;
  }

  int dim = 1 << n_modes;
  cmatrix_t *H = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    H->data[i] = c_zero();
  }

  for (int state = 0; state < dim; state++) {
    // Diagonal: on-site energies + nearest-neighbor interaction
    double diag = 0.0;

    for (int i = 0; i < n_modes; i++) {
      int bitpos = n_modes - 1 - i;
      int occ_i = (state >> bitpos) & 1;
      diag += epsilon[i] * occ_i;
    }

    for (int i = 0; i + 1 < n_modes; i++) {
      int occ_i = (state >> (n_modes - 1 - i)) & 1;
      int occ_ip1 = (state >> (n_modes - 1 - (i + 1))) & 1;
      diag += U * occ_i * occ_ip1;
    }

    CMAT(H, state, state) = c_add(CMAT(H, state, state), c_real(diag));

    // Hopping: -t * (a_i^\dagger a_{i+1} + a_{i+1}^\dagger a_i)
    for (int i = 0; i + 1 < n_modes; i++) {
      int sign1, sign2;

      // a_i^\dagger a_{i+1}: annihilate mode i+1, create mode i
      int mid = direct_annihilate(state, i + 1, n_modes, &sign1);
      if (mid >= 0) {
        int fin = direct_create(mid, i, n_modes, &sign2);

        if (fin >= 0) {
          double amp = -t * sign1 * sign2;
          CMAT(H, fin, state) = c_add(CMAT(H, fin, state), c_real(amp));
        }
      }

      // a_{i+1}^\dagger a_i: annihilate mode i, create mode i+1
      mid = direct_annihilate(state, i, n_modes, &sign1);
      if (mid >= 0) {
        int fin = direct_create(mid, i + 1, n_modes, &sign2);

        if (fin >= 0) {
          double amp = -t * sign1 * sign2;
          CMAT(H, fin, state) = c_add(CMAT(H, fin, state), c_real(amp));
        }
      }
    }
  }

  return H;
}
