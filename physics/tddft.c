#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "complex.h"
#include "dft.h"
#include "molecular_hf.h"
#include "physics/molecular_dft.h"
#include "physics/molecular_integrals.h"
#include "tddft.h"
#include <math.h>
#include <stdlib.h>

/*
 * Central finite-difference d(V_xc)/dn at density n, using the same
 * combined Slater+PZ81 potential molecular_ks_lda's SCF itself uses
 */
static double lda_fxc_kernel(double n) {
  if (n <= 0.0) {
    return 0.0;
  }

  double h = fmax(1e-6, n * 1e-5);
  double n_minus = n - h;
  if (n_minus < 0.0) {
    n_minus = 0.0;
  }

  return (lda_xc_potential(n + h) - lda_xc_potential(n_minus)) /
         (n + h - n_minus);
}

tddft_result_t *tddft_tda_lda(basis_function_t **basis, int n_basis,
                              const molecule_t *mol,
                              const molecular_dft_result_t *ks,
                              const molecular_grid_t *grid) {
  if (!basis || n_basis < 1 || !mol || !ks || !ks->converged || !grid ||
      grid->n_points <= 0 || ks->n_basis != n_basis) {
    return NULL;
  }

  if (ks->n_electrons < 0 || ks->n_electrons % 2 != 0 ||
      ks->n_electrons > 2 * n_basis) {
    return NULL;
  }

  int n_occ = ks->n_electrons / 2;
  int n_virt = n_basis - n_occ;
  if (n_occ <= 0 || n_virt <= 0) {
    return NULL;
  }

  int dim = n_occ * n_virt;
  int ng = grid->n_points;

  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, n_basis, mol);
  double *eri_ao = molecular_eri_tensor(basis, n_basis);
  double *h_mo = malloc((size_t)n_basis * n_basis * sizeof(double));
  double *eri_mo =
      malloc((size_t)n_basis * n_basis * n_basis * n_basis * sizeof(double));
  if (!h_ao || !eri_ao || !h_mo || !eri_mo) {
    cmatrix_free(h_ao);
    free(eri_ao);
    free(h_mo);
    free(eri_mo);

    return NULL;
  }

  molecular_ao_to_mo(h_ao, eri_ao, ks->C, n_basis, h_mo, eri_mo);
  cmatrix_free(h_ao);
  free(eri_ao);
  free(h_mo); /* TDA's A matrix only needs eri_mo; h_mo is unused here since
                 orbital-energy differences already fold in converged
                 core+Coulomb+xc contributions */

  /* MO values at every grid point: phi_p(r_g) = sum_mu C_mu_p chi_mu(r_g). */
  double *ao_vals = malloc((size_t)ng * n_basis * sizeof(double));
  double *mo_vals = malloc((size_t)ng * n_basis * sizeof(double));
  double *fxc_vals = malloc((size_t)ng * sizeof(double));
  if (!ao_vals || !mo_vals || !fxc_vals) {
    free(eri_mo);
    free(ao_vals);
    free(mo_vals);
    free(fxc_vals);

    return NULL;
  }

  for (int g = 0; g < ng; g++) {
    double r[3] = {grid->points[g].x, grid->points[g].y, grid->points[g].z};

    for (int p = 0; p < n_basis; p++) {
      ao_vals[g * n_basis + p] = basis_function_value(basis[p], r);
    }
  }

  for (int g = 0; g < ng; g++) {
    for (int p = 0; p < n_basis; p++) {
      double sum = 0.0;

      for (int mu = 0; mu < n_basis; mu++) {
        sum += ao_vals[g * n_basis + mu] * CMAT(ks->C, mu, p).re;
      }
      mo_vals[g * n_basis + p] = sum;
    }
  }

  free(ao_vals);

  /* Ground-state density (closed-shell: 2 electrons per occupied MO) and local
   * LDA xc kernel at each grid point */
  for (int g = 0; g < ng; g++) {
    double dens = 0.0;

    for (int i = 0; i < n_occ; i++) {
      // NOLINTNEXTLINE(clang-analyzer-core.uninitialized.Assign)
      double v = mo_vals[g * n_basis + i];

      dens += 2.0 * v * v;
    }
    fxc_vals[g] = lda_fxc_kernel(dens);
  }

  /* Build the TDA coupling matrix A (real, symmetric, dim x dim, indexed by
   * (i,a) -> i*n_virt + (a-n_occ)), then diagonalize via same complex Hermitian
   * solver used throughout this codebase for real symmetric matrices */
  cmatrix_t *A = cmatrix_alloc(dim, dim);
  if (!A) {
    free(eri_mo);
    free(mo_vals);
    free(fxc_vals);

    return NULL;
  }

  for (int i = 0; i < n_occ; i++) {
    for (int a = n_occ; a < n_basis; a++) {
      int row = i * n_virt + (a - n_occ);

      for (int j = 0; j < n_occ; j++) {
        for (int b = n_occ; b < n_basis; b++) {
          int col = j * n_virt + (b - n_occ);

          double val = 0.0;
          if (i == j && a == b) {
            val += ks->orbital_energies[a] - ks->orbital_energies[i];
          }
          val += 2.0 * MOLINT_ERI(eri_mo, n_basis, i, a, j, b);

          double fxc_int = 0.0;
          for (int g = 0; g < ng; g++) {
            fxc_int += grid->points[g].weight * fxc_vals[g] *
                       mo_vals[g * n_basis + i] * mo_vals[g * n_basis + a] *
                       mo_vals[g * n_basis + j] * mo_vals[g * n_basis + b];
          }
          val += 2.0 * fxc_int;

          CMAT(A, row, col) = c_real(val);
        }
      }
    }
  }

  free(eri_mo);
  free(mo_vals);
  free(fxc_vals);

  eigen_t *eig = cmatrix_eigh_complex(A);
  cmatrix_free(A);
  if (!eig) {
    return NULL;
  }

  tddft_result_t *res = malloc(sizeof(tddft_result_t));
  if (!res) {
    eigen_free(eig);

    return NULL;
  }

  res->excitation_energies = eig->eigenvalues; /* transfer ownership */
  res->dim = dim;
  res->n_occ = n_occ;
  res->n_virt = n_virt;
  res->n_basis = n_basis;

  cmatrix_free(eig->eigenvectors);
  free(eig);

  return res;
}

void tddft_result_free(tddft_result_t *res) {
  if (!res) {
    return;
  }

  free(res->excitation_energies);
  free(res);
}

