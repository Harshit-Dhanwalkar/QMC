#include "molecular_dft.h"
#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "dft.h"
#include "molecular_integrals.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------
 * Becke fuzzy-Voronoi molecular grid
 * ------------------------------------------------------------------- */

/*
 * Becke's iterated-polynomial smoothing step, applied k=3 times (Reference:
 * Becke 1988's own choice, standard in every subsequent implementation): each
 * application makes s(\mu) flatter near \mu=+-1 (i.e. sharper transition
 * region), giving faster convergence of the resulting cell functions.
 */
static double becke_smooth(double mu, int iterations) {
  double f = mu;
  for (int k = 0; k < iterations; k++) {
    f = 1.5 * f - 0.5 * f * f * f;
  }

  return f;
}

/*
 * Unnormalized cell weight P_i(r) for atom i (product over all j!=i of smoothed
 * step function s(\mu_ij(r))); `dist` is a caller-provided scratch array of
 * length mol->n_atoms holding |r - R_k| for every atom k (computed once per
 * grid point and reused across all i, avoiding n_atoms^3 repeated distance
 * computations).
 */
static double becke_cell_unnormalized(int i, const double *dist,
                                      const molecule_t *mol) {
  int n = mol->n_atoms;
  double P = 1.0;
  for (int j = 0; j < n; j++) {
    if (j == i) {
      continue;
    }

    double dx = mol->center[i][0] - mol->center[j][0];
    double dy = mol->center[i][1] - mol->center[j][1];
    double dz = mol->center[i][2] - mol->center[j][2];
    double Rij = sqrt(dx * dx + dy * dy + dz * dz);
    double mu = (dist[i] - dist[j]) / Rij;
    double f3 = becke_smooth(mu, 3);
    double s = 0.5 * (1.0 - f3);

    P *= s;
  }

  return P;
}

/*
 * Becke's normalized atomic partition weight w_i(r) = P_i(r) / \sum_k P_k(r).
 */
static double becke_weight_at(int i, const double r[3], const molecule_t *mol,
                              double *dist_scratch) {
  int n = mol->n_atoms;
  for (int k = 0; k < n; k++) {
    double dx = r[0] - mol->center[k][0];
    double dy = r[1] - mol->center[k][1];
    double dz = r[2] - mol->center[k][2];

    dist_scratch[k] = sqrt(dx * dx + dy * dy + dz * dz);
  }

  double Pi = becke_cell_unnormalized(i, dist_scratch, mol);
  if (Pi == 0.0) {
    return 0.0;
  }

  double sum = 0.0;
  for (int k = 0; k < n; k++) {
    sum += becke_cell_unnormalized(k, dist_scratch, mol);
  }

  return Pi / sum;
}

/*
 * Gauss-Chebyshev-of-the-second-kind radial grid, mapped (-1,1) -> (0, \infty)
 * via r = r_m*(1 + x)/(1 - x).
 * NOTE: The raw GC2 weight (\pi / (N + 1)) * \sin^2(\theta_i) carries an
 * implicit \sqrt(1 - x^2) = \sin(\theta_i) weight function, which must be
 * divided back out before multiplying by the r^2*dr/dx Jacobian to get a plain
 * (unweighted) integral over r in (0, \infty).
 */
static void gauss_chebyshev_radial(int n, double r_m, double *r, double *w) {
  for (int i = 1; i <= n; i++) {
    double theta = i * M_PI / (n + 1);
    double x = cos(theta);
    double ri = r_m * (1.0 + x) / (1.0 - x);
    double drdx = r_m * 2.0 / ((1.0 - x) * (1.0 - x));

    r[i - 1] = ri;
    w[i - 1] = (M_PI / (n + 1)) * sin(theta) * drdx * ri * ri;
  }
}

molecular_grid_t *molecular_grid_build(const molecule_t *mol, int n_radial,
                                       int n_polar, int n_azimuthal,
                                       double radial_scale) {
  if (!mol || mol->n_atoms <= 0 || n_radial <= 0 || n_polar <= 0 ||
      n_azimuthal <= 0 || radial_scale <= 0.0) {
    return NULL;
  }

  // Gauss-Legendre nodes/weights in cos(theta) in (-1,1)
  double *gl_x = malloc((size_t)n_polar * sizeof(double));
  double *gl_w = malloc((size_t)n_polar * sizeof(double));
  {
    // Newton's method on the Legendre polynomial
    for (int i = 0; i < n_polar; i++) {
      double x = cos(M_PI * (i + 0.75) / (n_polar + 0.5));
      double dpdx = 0.0;

      for (int iter = 0; iter < 100; iter++) {
        double p0 = 1.0, p1 = x;

        for (int k = 2; k <= n_polar; k++) {
          double p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
          p0 = p1;
          p1 = p2;
        }

        dpdx = n_polar * (x * p1 - p0) / (x * x - 1.0);
        double dx = -p1 / dpdx;

        x += dx;
        if (fabs(dx) < 1e-15) {
          break;
        }
      }

      gl_x[i] = x;
      gl_w[i] = 2.0 / ((1.0 - x * x) * dpdx * dpdx);
    }
  }

  double *r_nodes = malloc((size_t)n_radial * sizeof(double));
  double *r_weights = malloc((size_t)n_radial * sizeof(double));
  gauss_chebyshev_radial(n_radial, radial_scale, r_nodes, r_weights);

  int n_ang = n_polar * n_azimuthal;
  int points_per_atom = n_radial * n_ang;
  int total_points = mol->n_atoms * points_per_atom;

  molecular_grid_t *grid = malloc(sizeof(molecular_grid_t));
  grid->points = malloc((size_t)total_points * sizeof(dft_grid_point_t));
  grid->n_points = total_points;

  double *dist_scratch = malloc((size_t)mol->n_atoms * sizeof(double));
  double w_phi = 2.0 * M_PI / n_azimuthal;

  int idx = 0;
  for (int a = 0; a < mol->n_atoms; a++) {
    for (int kr = 0; kr < n_radial; kr++) {
      double r = r_nodes[kr];
      double wr = r_weights[kr];

      for (int kt = 0; kt < n_polar; kt++) {
        double ct = gl_x[kt];
        double st = sqrt(1.0 - ct * ct);
        double wt = gl_w[kt];

        for (int kp = 0; kp < n_azimuthal; kp++) {
          double phi = 2.0 * M_PI * kp / n_azimuthal;
          double pt[3] = {mol->center[a][0] + r * st * cos(phi),
                          mol->center[a][1] + r * st * sin(phi),
                          mol->center[a][2] + r * ct};

          double becke_w = becke_weight_at(a, pt, mol, dist_scratch);
          double combined = wr * wt * w_phi * becke_w;

          grid->points[idx].x = pt[0];
          grid->points[idx].y = pt[1];
          grid->points[idx].z = pt[2];
          grid->points[idx].weight = combined;

          idx++;
        }
      }
    }
  }

  free(dist_scratch);
  free(r_nodes);
  free(r_weights);
  free(gl_x);
  free(gl_w);

  return grid;
}

molecular_grid_t *molecular_grid_build_default(const molecule_t *mol) {
  return molecular_grid_build(mol, 60, 20, 30, 1.0);
}

void molecular_grid_free(molecular_grid_t *grid) {
  if (!grid) {
    return;
  }

  free(grid->points);
  free(grid);
}

/* ---------------------------------------------------------------------
 * Kohn-Sham LDA SCF
 * ------------------------------------------------------------------- */

molecular_dft_result_t *molecular_ks_lda(basis_function_t **basis, int n_basis,
                                         const molecule_t *mol, int n_electrons,
                                         const molecular_grid_t *grid,
                                         double mix, double conv_tol,
                                         int max_iter) {
  if (!basis || n_basis <= 0 || !mol || n_electrons <= 0 ||
      n_electrons % 2 != 0 || n_electrons > 2 * n_basis || !grid ||
      grid->n_points <= 0 || mix <= 0.0 || mix > 1.0 || max_iter < 1) {
    return NULL;
  }

  int n = n_basis;
  int n_occ = n_electrons / 2;
  int ng = grid->n_points;

  cmatrix_t *S = molecular_overlap_matrix(basis, n);
  cmatrix_t *Hcore = molecular_core_hamiltonian(basis, n, mol);
  double *eri = molecular_eri_tensor(basis, n);
  if (!S || !Hcore || !eri) {
    cmatrix_free(S);
    cmatrix_free(Hcore);
    free(eri);

    return NULL;
  }

  /* NOTE: Precompute every basis function's value at every grid point once
   * (dominant memory cost is n_basis*n_grid doubles: fine for small systems
   * this module targets; avoids re-evaluating primitives on every SCF
   * iteration). */
  double *ao_vals = malloc((size_t)ng * n * sizeof(double));
  for (int g = 0; g < ng; g++) {
    double r[3] = {grid->points[g].x, grid->points[g].y, grid->points[g].z};

    for (int p = 0; p < n; p++) {
      ao_vals[g * n + p] = basis_function_value(basis[p], r);
    }
  }

  // Lowdin symmetric orthogonalization
  eigen_t *eig_S = cmatrix_eigh_complex(S);
  cmatrix_t *Shalf = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;

      for (int k = 0; k < n; k++) {
        double inv_sqrt = 1.0 / sqrt(eig_S->eigenvalues[k]);

        sum += CMAT(eig_S->eigenvectors, i, k).re * inv_sqrt *
               CMAT(eig_S->eigenvectors, j, k).re;
      }

      CMAT(Shalf, i, j) = c_real(sum);
    }
  }

  eigen_free(eig_S);

  double *D = calloc((size_t)n * n, sizeof(double));
  double *D_computed = malloc((size_t)n * n * sizeof(double));
  double *D_new = malloc((size_t)n * n * sizeof(double));
  double *C_arr = malloc((size_t)n * n * sizeof(double));
  double *orbital_energies = malloc((size_t)n * sizeof(double));
  double *dens = malloc((size_t)ng * sizeof(double));
  double *vxc_pot = malloc((size_t)ng * sizeof(double));

  double E_total = 0.0, E_old = 0.0;
  double E_core = 0.0, E_coulomb = 0.0, E_xc = 0.0;
  int converged = 0;
  int iter = 0;

  for (; iter < max_iter; iter++) {
    // Classical (exact, analytic) Coulomb matrix J from the ERI tensor
    cmatrix_t *J = cmatrix_alloc(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double v = 0.0;

        for (int k = 0; k < n; k++) {
          for (int l = 0; l < n; l++) {
            v += D[k * n + l] * MOLINT_ERI(eri, n, i, j, l, k);
          }
        }

        CMAT(J, i, j) = c_real(v);
      }
    }

    // Density on the grid, from current (pre-update) density matrix
    for (int g = 0; g < ng; g++) {
      double v = 0.0;

      for (int p = 0; p < n; p++) {
        double ap = ao_vals[g * n + p];

        if (ap == 0.0) {
          continue;
        }

        for (int q = 0; q < n; q++) {
          v += D[p * n + q] * ap * ao_vals[g * n + q];
        }
      }

      dens[g] = v > 0.0 ? v : 0.0;
    }

    double Exc_g = 0.0, Vxc_n_g = 0.0;
    for (int g = 0; g < ng; g++) {
      double n_g = dens[g] > 1e-12 ? dens[g] : 1e-12;
      double eps_xc = lda_xc_energy_density(n_g);

      vxc_pot[g] = lda_xc_potential(n_g);
      /* Weight the energy/potential contributions by unfloored density, so an
       * exactly-zero-density point contributes exactly zero rather than a
       * phantom floor-value times a (possibly large, far-radial-shell)
       * quadrature weight. */
      Exc_g += grid->points[g].weight * dens[g] * eps_xc;
      Vxc_n_g += grid->points[g].weight * dens[g] * vxc_pot[g];
    }

    // Vxc matrix, built on the grid
    cmatrix_t *Vxc = cmatrix_alloc(n, n);
    for (int p = 0; p < n; p++) {
      for (int q = p; q < n; q++) {
        double v = 0.0;

        for (int g = 0; g < ng; g++) {
          v += grid->points[g].weight * vxc_pot[g] * ao_vals[g * n + p] *
               ao_vals[g * n + q];
        }

        CMAT(Vxc, p, q) = c_real(v);
        CMAT(Vxc, q, p) = c_real(v);
      }
    }

    // F = Hcore + J + Vxc, Lowdin-orthogonalize, diagonalize
    cmatrix_t *F = cmatrix_alloc(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        CMAT(F, i, j) = c_real(CMAT(Hcore, i, j).re + CMAT(J, i, j).re +
                               CMAT(Vxc, i, j).re);
      }
    }

    cmatrix_t *tmp = cmatrix_multiply(Shalf, F);
    cmatrix_t *Fp = cmatrix_multiply(tmp, Shalf);

    cmatrix_free(tmp);
    cmatrix_free(F);

    eigen_t *eig_F = cmatrix_eigh_complex(Fp);
    cmatrix_free(Fp);

    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        double v = 0.0;

        for (int p = 0; p < n; p++) {
          v += CMAT(Shalf, i, p).re * CMAT(eig_F->eigenvectors, p, k).re;
        }

        C_arr[i * n + k] = v;
      }
    }

    memcpy(orbital_energies, eig_F->eigenvalues, (size_t)n * sizeof(double));

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double v = 0.0;

        for (int k = 0; k < n_occ; k++) {
          v += 2.0 * C_arr[i * n + k] * C_arr[j * n + k];
        }

        D_computed[i * n + j] = v;
      }
    }

    /* NOTE: Linear density mixing/damping: plain unmixed SCF converges fine for
     * simplest case (H2) but oscillates for anything with near-degenerate
     * orbitals (e.g. LiH) */
    for (int i = 0; i < n * n; i++) {
      D_new[i] = mix * D_computed[i] + (1.0 - mix) * D[i];
    }

    double E_J = 0.0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        E_J += D_new[i * n + j] * CMAT(J, i, j).re;
      }
    }

    double sum_eps = 0.0;
    for (int i = 0; i < n_occ; i++) {
      sum_eps += eig_F->eigenvalues[i];
    }

    // KS total-energy formula
    double E_elec = 2.0 * sum_eps - 0.5 * E_J - Vxc_n_g + Exc_g;
    double E_nuc = molecule_nuclear_repulsion(mol);

    E_total = E_elec + E_nuc;
    E_core = 0.0;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        E_core += D_new[i * n + j] * CMAT(Hcore, i, j).re;
      }
    }

    E_coulomb = 0.5 * E_J;
    E_xc = Exc_g;

    eigen_free(eig_F);
    cmatrix_free(J);
    cmatrix_free(Vxc);

    memcpy(D, D_new, (size_t)n * n * sizeof(double));

    if (iter > 5 && fabs(E_total - E_old) < conv_tol) {
      converged = 1;
      iter++;

      break;
    }

    E_old = E_total;
  }

  molecular_dft_result_t *res = malloc(sizeof(molecular_dft_result_t));
  res->total_energy = E_total;
  res->e_core = E_core;
  res->e_coulomb = E_coulomb;
  res->e_xc = E_xc;
  res->e_nuclear = molecule_nuclear_repulsion(mol);
  res->orbital_energies = orbital_energies;
  res->n_basis = n;
  res->n_electrons = n_electrons;
  res->converged = converged;
  res->iterations = iter;

  res->C = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(res->C, i, j) = c_real(C_arr[i * n + j]);
    }
  }

  free(D);
  free(D_computed);
  free(D_new);
  free(C_arr);
  free(dens);
  free(vxc_pot);
  free(ao_vals);
  cmatrix_free(S);
  cmatrix_free(Hcore);
  cmatrix_free(Shalf);
  free(eri);

  return res;
}

molecular_dft_result_t *molecular_ks_lda_default(basis_function_t **basis,
                                                 int n_basis,
                                                 const molecule_t *mol,
                                                 int n_electrons,
                                                 const molecular_grid_t *grid) {
  return molecular_ks_lda(basis, n_basis, mol, n_electrons, grid, 0.3, 1e-9,
                          200);
}

/* ---------------------------------------------------------------------
 * Kohn-Sham PBE (GGA) SCF
 * ------------------------------------------------------------------- */
molecular_dft_result_t *molecular_ks_pbe(basis_function_t **basis, int n_basis,
                                         const molecule_t *mol, int n_electrons,
                                         const molecular_grid_t *grid,
                                         double mix, double conv_tol,
                                         int max_iter) {
  if (!basis || n_basis <= 0 || !mol || n_electrons <= 0 ||
      n_electrons % 2 != 0 || n_electrons > 2 * n_basis || !grid ||
      grid->n_points <= 0 || mix <= 0.0 || mix > 1.0 || max_iter < 1) {
    return NULL;
  }

  int n = n_basis;
  int n_occ = n_electrons / 2;
  int ng = grid->n_points;

  cmatrix_t *S = molecular_overlap_matrix(basis, n);
  cmatrix_t *Hcore = molecular_core_hamiltonian(basis, n, mol);
  double *eri = molecular_eri_tensor(basis, n);
  if (!S || !Hcore || !eri) {
    cmatrix_free(S);
    cmatrix_free(Hcore);
    free(eri);

    return NULL;
  }

  /* Precompute every basis function's value and gradient at every grid point
   * once.
   */
  double *ao_vals = malloc((size_t)ng * n * sizeof(double));
  double *ao_grads = malloc((size_t)ng * n * 3 * sizeof(double));
  for (int g = 0; g < ng; g++) {
    double r[3] = {grid->points[g].x, grid->points[g].y, grid->points[g].z};

    for (int p = 0; p < n; p++) {
      ao_vals[g * n + p] = basis_function_value(basis[p], r);
      basis_function_gradient(basis[p], r, &ao_grads[(g * n + p) * 3]);
    }
  }

  // Lowdin symmetric orthogonalization
  eigen_t *eig_S = cmatrix_eigh_complex(S);
  cmatrix_t *Shalf = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;

      for (int k = 0; k < n; k++) {
        double inv_sqrt = 1.0 / sqrt(eig_S->eigenvalues[k]);

        sum += CMAT(eig_S->eigenvectors, i, k).re * inv_sqrt *
               CMAT(eig_S->eigenvectors, j, k).re;
      }

      CMAT(Shalf, i, j) = c_real(sum);
    }
  }

  eigen_free(eig_S);

  double *D = calloc((size_t)n * n, sizeof(double));
  double *D_computed = malloc((size_t)n * n * sizeof(double));
  double *D_new = malloc((size_t)n * n * sizeof(double));
  double *C_arr = malloc((size_t)n * n * sizeof(double));
  double *orbital_energies = malloc((size_t)n * sizeof(double));
  double *dens = malloc((size_t)ng * sizeof(double));
  double *dens_grad = malloc((size_t)ng * 3 * sizeof(double));
  double *vrho_pot = malloc((size_t)ng * sizeof(double));
  double *vsigma_pot = malloc((size_t)ng * sizeof(double));

  double E_total = 0.0, E_old = 0.0;
  double E_core = 0.0, E_coulomb = 0.0, E_xc = 0.0;
  int converged = 0;
  int iter = 0;

  for (; iter < max_iter; iter++) {
    // Classical (exact, analytic) Coulomb matrix J from ERI tensor
    cmatrix_t *J = cmatrix_alloc(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double v = 0.0;

        for (int k = 0; k < n; k++) {
          for (int l = 0; l < n; l++) {
            v += D[k * n + l] * MOLINT_ERI(eri, n, i, j, l, k);
          }
        }

        CMAT(J, i, j) = c_real(v);
      }
    }

    // Density AND density gradient on the grid, from the current (pre-update)
    // density matrix:
    // n(r) = \sum_{pq} D_{pq} * \phi_p * \phi_q, grad n(r)
    //      = 2 * \sum_{pq} D_{pq} * \grad(\phi_p) * \phi_q
    // (D symmetric)
    for (int g = 0; g < ng; g++) {
      double v = 0.0;
      double gx = 0.0, gy = 0.0, gz = 0.0;

      for (int p = 0; p < n; p++) {
        double ap = ao_vals[g * n + p];
        const double *gp = &ao_grads[(g * n + p) * 3];

        for (int q = 0; q < n; q++) {
          double Dpq = D[p * n + q];
          double aq = ao_vals[g * n + q];

          v += Dpq * ap * aq;
          gx += Dpq * gp[0] * aq;
          gy += Dpq * gp[1] * aq;
          gz += Dpq * gp[2] * aq;
        }
      }

      dens[g] = v > 0.0 ? v : 0.0;
      dens_grad[g * 3 + 0] = 2.0 * gx;
      dens_grad[g * 3 + 1] = 2.0 * gy;
      dens_grad[g * 3 + 2] = 2.0 * gz;
    }

    double Exc_g = 0.0;
    for (int g = 0; g < ng; g++) {
      double n_g = dens[g] > 1e-12 ? dens[g] : 1e-12;
      double dgx = dens_grad[g * 3 + 0], dgy = dens_grad[g * 3 + 1],
             dgz = dens_grad[g * 3 + 2];
      double sigma_g = dgx * dgx + dgy * dgy + dgz * dgz;

      double eps_xc = pbe_xc_energy_density(n_g, sigma_g);
      double vrho;
      double vsigma;

      pbe_xc_potential(n_g, sigma_g, &vrho, &vsigma);

      vrho_pot[g] = vrho;
      vsigma_pot[g] = vsigma;
      // Weight by unfloored density
      Exc_g += grid->points[g].weight * dens[g] * eps_xc;
    }

    // GGA Vxc matrix: V_pq = \int [ vrho * \phi_p * \phi_q + 2 * vsigma *
    // (\grad n).(\phi_q * \grad(\phi_p) + \phi_p * \grad(\phi_q)) ] dr
    cmatrix_t *Vxc = cmatrix_alloc(n, n);
    for (int p = 0; p < n; p++) {
      for (int q = p; q < n; q++) {
        double v = 0.0;

        for (int g = 0; g < ng; g++) {
          double ap = ao_vals[g * n + p], aq = ao_vals[g * n + q];
          const double *gp = &ao_grads[(g * n + p) * 3];
          const double *gq = &ao_grads[(g * n + q) * 3];
          double dgx = dens_grad[g * 3 + 0], dgy = dens_grad[g * 3 + 1],
                 dgz = dens_grad[g * 3 + 2];

          double grad_dot = dgx * (aq * gp[0] + ap * gq[0]) +
                            dgy * (aq * gp[1] + ap * gq[1]) +
                            dgz * (aq * gp[2] + ap * gq[2]);

          v += grid->points[g].weight *
               (vrho_pot[g] * ap * aq + 2.0 * vsigma_pot[g] * grad_dot);
        }

        CMAT(Vxc, p, q) = c_real(v);
        CMAT(Vxc, q, p) = c_real(v);
      }
    }

    // F = Hcore + J + Vxc, Lowdin-orthogonalize, diagonalize
    cmatrix_t *F = cmatrix_alloc(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        CMAT(F, i, j) = c_real(CMAT(Hcore, i, j).re + CMAT(J, i, j).re +
                               CMAT(Vxc, i, j).re);
      }
    }

    cmatrix_t *tmp = cmatrix_multiply(Shalf, F);
    cmatrix_t *Fp = cmatrix_multiply(tmp, Shalf);

    cmatrix_free(tmp);
    cmatrix_free(F);

    eigen_t *eig_F = cmatrix_eigh_complex(Fp);
    cmatrix_free(Fp);

    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        double v = 0.0;

        for (int p = 0; p < n; p++) {
          v += CMAT(Shalf, i, p).re * CMAT(eig_F->eigenvectors, p, k).re;
        }

        C_arr[i * n + k] = v;
      }
    }

    memcpy(orbital_energies, eig_F->eigenvalues, (size_t)n * sizeof(double));

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double v = 0.0;

        for (int k = 0; k < n_occ; k++) {
          v += 2.0 * C_arr[i * n + k] * C_arr[j * n + k];
        }

        D_computed[i * n + j] = v;
      }
    }

    // Linear density mixing
    for (int i = 0; i < n * n; i++) {
      D_new[i] = mix * D_computed[i] + (1.0 - mix) * D[i];
    }

    double E_J = 0.0, Vxc_n_g = 0.0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        E_J += D_new[i * n + j] * CMAT(J, i, j).re;
        Vxc_n_g += D_new[i * n + j] * CMAT(Vxc, i, j).re;
      }
    }

    double sum_eps = 0.0;
    for (int i = 0; i < n_occ; i++) {
      sum_eps += eig_F->eigenvalues[i];
    }

    // KS total-energy formula
    double E_elec = 2.0 * sum_eps - 0.5 * E_J - Vxc_n_g + Exc_g;
    double E_nuc = molecule_nuclear_repulsion(mol);

    E_total = E_elec + E_nuc;
    E_core = 0.0;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        E_core += D_new[i * n + j] * CMAT(Hcore, i, j).re;
      }
    }

    E_coulomb = 0.5 * E_J;
    E_xc = Exc_g;

    eigen_free(eig_F);
    cmatrix_free(J);
    cmatrix_free(Vxc);

    memcpy(D, D_new, (size_t)n * n * sizeof(double));

    if (iter > 5 && fabs(E_total - E_old) < conv_tol) {
      converged = 1;
      iter++;

      break;
    }

    E_old = E_total;
  }

  molecular_dft_result_t *res = malloc(sizeof(molecular_dft_result_t));

  res->total_energy = E_total;
  res->e_core = E_core;
  res->e_coulomb = E_coulomb;
  res->e_xc = E_xc;
  res->e_nuclear = molecule_nuclear_repulsion(mol);
  res->orbital_energies = orbital_energies;
  res->n_basis = n;
  res->n_electrons = n_electrons;
  res->converged = converged;
  res->iterations = iter;

  res->C = cmatrix_alloc(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(res->C, i, j) = c_real(C_arr[i * n + j]);
    }
  }

  free(D);
  free(D_computed);
  free(D_new);
  free(C_arr);
  free(dens);
  free(dens_grad);
  free(vrho_pot);
  free(vsigma_pot);
  free(ao_vals);
  free(ao_grads);
  cmatrix_free(S);
  cmatrix_free(Hcore);
  cmatrix_free(Shalf);
  free(eri);

  return res;
}

molecular_dft_result_t *molecular_ks_pbe_default(basis_function_t **basis,
                                                 int n_basis,
                                                 const molecule_t *mol,
                                                 int n_electrons,
                                                 const molecular_grid_t *grid) {
  return molecular_ks_pbe(basis, n_basis, mol, n_electrons, grid, 0.3, 1e-9,
                          200);
}

void molecular_dft_result_free(molecular_dft_result_t *res) {
  if (!res) {
    return;
  }

  free(res->orbital_energies);
  cmatrix_free(res->C);
  free(res);
}
