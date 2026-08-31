#include "fci.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "second_quant.h"
#include <stdlib.h>

static int binomial(int n, int k) {
  if (k < 0 || k > n) {
    return 0;
  }
  if (k > n - k) {
    k = n - k;
  }

  long long result = 1;
  for (int i = 0; i < k; i++) {
    result = result * (n - i) / (i + 1);
  }

  return (int)result;
}

/* Shared implementation: takes an already-built full Fock-space Hamiltonian
 * (dimension full_dim = 2^n_modes) and extracts/diagonalizes the
 * n_electrons-particle sector
 */
static fci_result_t *fci_solve_from_full_hamiltonian(cmatrix_t *H, int n_modes,
                                                     int n_electrons,
                                                     int n_spatial_report) {
  if (!H) {
    return NULL;
  }

  int full_dim = H->nrows;
  if (n_electrons < 0 || n_electrons > n_modes) {
    cmatrix_free(H);

    return NULL;
  }

  int dim = binomial(n_modes, n_electrons);
  int *basis_states = malloc((size_t)dim * sizeof(int));
  if (!basis_states) {
    cmatrix_free(H);

    return NULL;
  }

  int count = 0;
  for (int state = 0; state < full_dim; state++) {
    if (__builtin_popcount((unsigned)state) == n_electrons) {
      basis_states[count++] = state;
    }
  }

  /* count must equal dim by combinatorics of popcount; if it doesn't, something
   * is inconsistent (e.g. n_modes not matching H's actual mode count) - fail
   * cleanly. */
  if (count != dim) {
    free(basis_states);
    cmatrix_free(H);

    return NULL;
  }

  cmatrix_t *H_sector = cmatrix_alloc(dim, dim);
  if (!H_sector) {
    free(basis_states);
    cmatrix_free(H);

    return NULL;
  }

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      CMAT(H_sector, i, j) = CMAT(H, basis_states[i], basis_states[j]);
    }
  }

  cmatrix_free(H);

  eigen_t *eig = cmatrix_eigh_complex(H_sector);
  cmatrix_free(H_sector);
  if (!eig) {
    free(basis_states);

    return NULL;
  }

  fci_result_t *res = malloc(sizeof(fci_result_t));
  if (!res) {
    free(basis_states);
    eigen_free(eig);

    return NULL;
  }

  res->dim = dim;
  res->n_spatial = n_spatial_report;
  res->n_electrons = n_electrons;
  res->basis_states = basis_states;
  res->eigenvalues = eig->eigenvalues;   /* transfer ownership */
  res->eigenvectors = eig->eigenvectors; /* transfer ownership */
  res->ground_energy = dim > 0 ? eig->eigenvalues[0] : 0.0;

  free(eig); // frees only eigen_t wrapper, not arrays

  return res;
}

fci_result_t *fci_solve(int n_spatial, const double *h_mo, const double *eri_mo,
                        double nuclear_repulsion, int n_electrons) {
  if (n_spatial < 1 || !h_mo || !eri_mo) {
    return NULL;
  }

  int n_modes = 2 * n_spatial;
  if (n_electrons < 0 || n_electrons > n_modes) {
    return NULL;
  }

  cmatrix_t *H = second_quant_build_molecular_hamiltonian(
      n_spatial, h_mo, eri_mo, nuclear_repulsion);

  return fci_solve_from_full_hamiltonian(H, n_modes, n_electrons, n_spatial);
}

fci_result_t *fci_solve_frozen_core(int n_spatial, int n_frozen,
                                    const double *h_mo, const double *eri_mo,
                                    double nuclear_repulsion,
                                    int n_electrons_active) {
  if (n_spatial < 1 || n_frozen < 0 || n_frozen >= n_spatial || !h_mo ||
      !eri_mo) {
    return NULL;
  }

  int n_active_spatial = n_spatial - n_frozen;
  int n_modes = 2 * n_active_spatial;
  if (n_electrons_active < 0 || n_electrons_active > n_modes) {
    return NULL;
  }

  cmatrix_t *H = second_quant_build_molecular_hamiltonian_frozen_core(
      n_spatial, n_frozen, h_mo, eri_mo, nuclear_repulsion);
  return fci_solve_from_full_hamiltonian(H, n_modes, n_electrons_active,
                                         n_active_spatial);
}

void fci_result_free(fci_result_t *res) {
  if (!res) {
    return;
  }

  free(res->eigenvalues);
  cmatrix_free(res->eigenvectors);
  free(res->basis_states);
  free(res);
}
