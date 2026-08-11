/*
Restricted Hartree-Fock for closed-shell, s-orbitals-only atoms/ions.
*/

#include "hartree_fock.h"
#include "../core/complex.h"
#include "../core/linalg/linalg.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "central_potential.h"
#include "potentials.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Rectangle-rule quadrature normalization: \int u^2 dr = 1.
static void normalize_u(double *u, int N, double dr) {
  double norm_sq = 0.0;
  for (int i = 0; i < N; i++) {
    norm_sq += u[i] * u[i] * dr;
  }

  if (norm_sq < 1e-300) {
    return;
  }

  double norm = sqrt(norm_sq);
  for (int i = 0; i < N; i++) {
    u[i] /= norm;
  }
}

/*
 * Y0_{ab}(r_i) = (1/r_i) * \int_0^{r_i} u_a(r')u_b(r') dr' +
 *                \int_{r_i}^{r_max} [u_a(r')u_b(r')/r'] dr'
 * NOTE: The l=0 ("monopole") radial Coulomb kernel: electrostatic potential of
 * a spherical charge density proportional to u_a(r)u_b(r), and (for a=b)
 * reduces to "charge enclosed / r + charge outside" shell-theorem formula.
 */
void compute_Y0(const double *r, int N, double dr, const double *ua,
                const double *ub, double *Y0_out) {
  double *inner = malloc(N * sizeof *inner);
  double *outer = malloc(N * sizeof *outer);
  if (!inner || !outer) {
    free(inner);
    free(outer);

    for (int i = 0; i < N; i++) {
      Y0_out[i] = 0.0;
    }

    return;
  }

  inner[0] = 0.0;
  for (int i = 1; i < N; i++) {
    double f_prev = ua[i - 1] * ub[i - 1];
    double f_curr = ua[i] * ub[i];

    inner[i] = inner[i - 1] + 0.5 * (f_prev + f_curr) * dr;
  }

  outer[N - 1] = 0.0;
  for (int i = N - 2; i >= 0; i--) {
    double f_prev = (r[i] > 0.0) ? ua[i] * ub[i] / r[i] : 0.0;
    double f_curr = (r[i + 1] > 0.0) ? ua[i + 1] * ub[i + 1] / r[i + 1] : 0.0;

    outer[i] = outer[i + 1] + 0.5 * (f_prev + f_curr) * dr;
  }

  for (int i = 0; i < N; i++) {
    double inv_r = (r[i] > 0.0) ? 1.0 / r[i] : 0.0;
    Y0_out[i] = inv_r * inner[i] + outer[i];
  }

  free(inner);
  free(outer);
}

/*
 * Dense real-symmetric (stored as complex-Hermitian, imaginary parts 0)
 * Fock matrix on the radial grid:
 *   F = T + V_nuc + \sum_k 2 * Y0_kk(r) * \delta_ab - \sum_k K_k(a,b)
 *
 *  Where
 *   K_k(a,b) = dr * u_k(r_a) * u_k(r_b) / max(r_a, r_b)
 * with Dirichlet boundary conditions u(r[0]) = u(r[N-1]) = 0 enforced by fully
 * decoupling boundary rows/columns.
 */
static cmatrix_t *build_fock_matrix(const double *r, int N, double dr, double Z,
                                    double **u, int n_orbitals) {
  double coeff = 0.5; // \hbar^2/(2m), atomic units: \hbar=m=1
  double diag_factor = 2.0 * coeff / (dr * dr);
  double offdiag_factor = -coeff / (dr * dr);

  cmatrix_t *F = cmatrix_alloc(N, N);
  if (!F) {
    return NULL;
  }

  for (int i = 0; i < N; i++) {
    double Vnuc = (r[i] > 0.0) ? -Z / r[i] : 0.0;
    CMAT(F, i, i) = c_real(diag_factor + Vnuc);

    if (i > 0) {
      CMAT(F, i, i - 1) = c_real(offdiag_factor);
    }
    if (i < N - 1) {
      CMAT(F, i, i + 1) = c_real(offdiag_factor);
    }
  }

  double *Y0 = malloc(N * sizeof *Y0);
  if (!Y0) {
    cmatrix_free(F);

    return NULL;
  }

  for (int k = 0; k < n_orbitals; k++) {
    compute_Y0(r, N, dr, u[k], u[k], Y0);

    for (int i = 0; i < N; i++) {
      CMAT(F, i, i) = c_add(CMAT(F, i, i), c_real(2.0 * Y0[i]));
    }
  }
  free(Y0);

  for (int k = 0; k < n_orbitals; k++) {
    for (int a = 0; a < N; a++) {
      for (int b = 0; b < N; b++) {
        double rmax_ab = (r[a] > r[b]) ? r[a] : r[b];
        if (rmax_ab <= 0.0) {
          continue;
        }

        double Kab = dr * u[k][a] * u[k][b] / rmax_ab;
        CMAT(F, a, b) = c_sub(CMAT(F, a, b), c_real(Kab));
      }
    }
  }

  // Dirichlet boundaries
  // NOTE: fully decouple endpoints so u(r_min)=u(r_max)=0 exactly, regardless
  // of dense exchange contribution above.
  for (int j = 0; j < N; j++) {
    CMAT(F, 0, j) = c_zero();
    CMAT(F, N - 1, j) = c_zero();
    CMAT(F, j, 0) = c_zero();
    CMAT(F, j, N - 1) = c_zero();
  }

  CMAT(F, 0, 0) = c_real(1e6 * diag_factor);
  CMAT(F, N - 1, N - 1) = c_real(1e6 * diag_factor);

  return F;
}

static void free_orbital_arrays(double **u, int n_orbitals) {
  if (!u) {
    return;
  }

  for (int k = 0; k < n_orbitals; k++) {
    free(u[k]);
  }

  free(u);
}

hf_result_t *hartree_fock_atom_s_orbitals(double *r, int N, double Z,
                                          int n_orbitals, double mix,
                                          double tol, int max_iter) {
  if (!r || N < 10 || Z <= 0.0 || n_orbitals < 1 || mix <= 0.0 || mix > 1.0 ||
      max_iter < 1) {
    return NULL;
  }

  double dr = r[1] - r[0];
  if (dr <= 0.0) {
    return NULL;
  }

  double **u = malloc(n_orbitals * sizeof(double *));
  if (!u) {
    return NULL;
  }

  for (int k = 0; k < n_orbitals; k++) {
    u[k] = malloc(N * sizeof(double));
  }

  // Initial guess: lowest n_orbitals bare-nucleus (l=0) eigenstates
  double Zc = Z;
  eigen_t *eig0 =
      central_potential_radial_solve(r, N, 0, 1.0, 1.0, V_coulomb, &Zc);

  if (!eig0 || eig0->n < n_orbitals) {
    if (eig0) {
      eigen_free(eig0);
    }

    free_orbital_arrays(u, n_orbitals);

    return NULL;
  }

  for (int k = 0; k < n_orbitals; k++) {
    for (int i = 0; i < N; i++) {
      u[k][i] = CMAT(eig0->eigenvectors, i, k).re;
    }

    normalize_u(u[k], N, dr);
  }

  eigen_free(eig0);

  double *orbital_energies = malloc(n_orbitals * sizeof(double));
  if (!orbital_energies) {
    free_orbital_arrays(u, n_orbitals);

    return NULL;
  }

  int converged = 0;
  int iter = 0;
  for (; iter < max_iter; iter++) {
    cmatrix_t *F = build_fock_matrix(r, N, dr, Z, u, n_orbitals);

    if (!F) {
      free_orbital_arrays(u, n_orbitals);
      free(orbital_energies);

      return NULL;
    }

    eigen_t *eig = cmatrix_eigh_generic(F);
    cmatrix_free(F);
    if (!eig || eig->n < n_orbitals) {
      if (eig) {
        eigen_free(eig);
      }

      free_orbital_arrays(u, n_orbitals);
      free(orbital_energies);

      return NULL;
    }

    double max_delta = 0.0;
    for (int k = 0; k < n_orbitals; k++) {
      double *u_new = malloc(N * sizeof(double));
      for (int i = 0; i < N; i++) {
        u_new[i] = CMAT(eig->eigenvectors, i, k).re;
      }

      normalize_u(u_new, N, dr);

      // Fix arbitrary overall sign of eigenvector before mixing.
      double dot = 0.0;
      for (int i = 0; i < N; i++) {
        dot += u[k][i] * u_new[i] * dr;
      }

      if (dot < 0.0) {
        for (int i = 0; i < N; i++) {
          u_new[i] = -u_new[i];
        }
      }

      for (int i = 0; i < N; i++) {
        double mixed = (1.0 - mix) * u[k][i] + mix * u_new[i];
        double delta = fabs(mixed - u[k][i]);

        if (delta > max_delta) {
          max_delta = delta;
        }

        u[k][i] = mixed;
      }

      normalize_u(u[k], N, dr);
      orbital_energies[k] = eig->eigenvalues[k];

      free(u_new);
    }

    eigen_free(eig);

    if (max_delta < tol) {
      converged = 1;
      iter++;
      break;
    }
  }

  // Fock and eigenpairs pass on converged orbitals
  // NOTE: so that returned energies/orbitals are mutually an eigenpair of
  // same Fock matrix.
  cmatrix_t *F_final = build_fock_matrix(r, N, dr, Z, u, n_orbitals);
  eigen_t *eig_final = F_final ? cmatrix_eigh_generic(F_final) : NULL;
  if (F_final) {
    cmatrix_free(F_final);
  }

  if (!eig_final || eig_final->n < n_orbitals) {
    if (eig_final) {
      eigen_free(eig_final);
    }

    free_orbital_arrays(u, n_orbitals);
    free(orbital_energies);

    return NULL;
  }

  for (int k = 0; k < n_orbitals; k++) {
    for (int i = 0; i < N; i++) {
      u[k][i] = CMAT(eig_final->eigenvectors, i, k).re;
    }

    normalize_u(u[k], N, dr);
    orbital_energies[k] = eig_final->eigenvalues[k];
  }

  int n_virtual = N - n_orbitals;
  double *virtual_energies =
      malloc((size_t)n_virtual * sizeof *virtual_energies);
  cvector_t **virtual_orbitals =
      malloc((size_t)n_virtual * sizeof *virtual_orbitals);

  if (!virtual_energies || !virtual_orbitals) {
    free(virtual_energies);
    free(virtual_orbitals);
    eigen_free(eig_final);
    free_orbital_arrays(u, n_orbitals);
    free(orbital_energies);

    return NULL;
  }

  for (int k = 0; k < n_virtual; k++) {
    int col = n_orbitals + k;
    double *uv = malloc((size_t)N * sizeof *uv);
    for (int i = 0; i < N; i++) {
      uv[i] = CMAT(eig_final->eigenvectors, i, col).re;
    }

    normalize_u(uv, N, dr);

    virtual_orbitals[k] = cvector_alloc(N);
    for (int i = 0; i < N; i++) {
      virtual_orbitals[k]->data[i] = c_real(uv[i]);
    }

    free(uv);
    virtual_energies[k] = eig_final->eigenvalues[col];
  }

  eigen_free(eig_final);

  // Total electronic energy: E = \sum_k 2 * eps_k - \sum_{i,j} (2*J_ij - K_ij)
  double total_energy = 0.0;
  for (int k = 0; k < n_orbitals; k++) {
    total_energy += 2.0 * orbital_energies[k];
  }

  double *Y0_buf = malloc(N * sizeof *Y0_buf);
  double correction = 0.0;
  for (int i = 0; i < n_orbitals; i++) {
    for (int j = 0; j < n_orbitals; j++) {
      compute_Y0(r, N, dr, u[j], u[j], Y0_buf);
      double Jij = 0.0;

      for (int a = 0; a < N; a++) {
        Jij += u[i][a] * u[i][a] * Y0_buf[a] * dr;
      }

      compute_Y0(r, N, dr, u[i], u[j], Y0_buf);
      double Kij = 0.0;
      for (int a = 0; a < N; a++) {
        Kij += u[i][a] * u[j][a] * Y0_buf[a] * dr;
      }

      correction += 2.0 * Jij - Kij;
    }
  }

  free(Y0_buf);
  total_energy -= correction;

  hf_result_t *res = malloc(sizeof(hf_result_t));
  if (!res) {
    free_orbital_arrays(u, n_orbitals);
    free(orbital_energies);

    return NULL;
  }

  res->n_orbitals = n_orbitals;
  res->N = N;
  res->Z = Z;
  res->orbital_energies = orbital_energies;
  res->total_energy = total_energy;
  res->iterations = iter;
  res->converged = converged;
  res->orbitals = malloc(n_orbitals * sizeof(cvector_t *));

  for (int k = 0; k < n_orbitals; k++) {
    res->orbitals[k] = cvector_alloc(N);

    for (int i = 0; i < N; i++) {
      res->orbitals[k]->data[i] = c_real(u[k][i]);
    }
  }
  free_orbital_arrays(u, n_orbitals);

  res->n_virtual = n_virtual;
  res->virtual_energies = virtual_energies;
  res->virtual_orbitals = virtual_orbitals;

  return res;
}

void hf_result_free(hf_result_t *res) {
  if (!res) {
    return;
  }

  free(res->orbital_energies);
  if (res->orbitals) {
    for (int k = 0; k < res->n_orbitals; k++) {
      cvector_free(res->orbitals[k]);
    }

    free(res->orbitals);
  }

  free(res->virtual_energies);
  if (res->virtual_orbitals) {
    for (int k = 0; k < res->n_virtual; k++) {
      cvector_free(res->virtual_orbitals[k]);
    }

    free(res->virtual_orbitals);
  }

  free(res);
}
