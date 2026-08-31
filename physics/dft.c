#include "dft.h"
#include "../core/complex.h"
#include "../core/linalg/eigen_generic.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "central_potential.h"
#include "hartree_fock.h" // reuse compute_Y0: same l=0 Coulomb kernel
#include "potentials.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------
 * LDA exchange-correlation functional
 * ------------------------------------------------------------------- */

double lda_exchange_energy_density(double n) {
  if (n <= 0.0) {
    return 0.0;
  }

  return -(3.0 / 4.0) * pow(3.0 / M_PI, 1.0 / 3.0) * pow(n, 1.0 / 3.0);
}

double lda_exchange_potential(double n) {
  if (n <= 0.0) {
    return 0.0;
  }

  return -pow(3.0 / M_PI, 1.0 / 3.0) * pow(n, 1.0 / 3.0);
}

static double rs_of_density(double n) {
  return pow(3.0 / (4.0 * M_PI * n), 1.0 / 3.0);
}

double lda_correlation_energy_density_pz81(double n) {
  if (n <= 0.0) {
    return 0.0;
  }

  double rs = rs_of_density(n);
  if (rs < 1.0) {
    const double A = 0.0311, B = -0.048, C = 0.0020, D = -0.0116;

    return A * log(rs) + B + C * rs * log(rs) + D * rs;
  } else {
    const double gamma = -0.1423, beta1 = 1.0529, beta2 = 0.3334;

    return gamma / (1.0 + beta1 * sqrt(rs) + beta2 * rs);
  }
}

double lda_correlation_potential_pz81(double n) {
  if (n <= 0.0) {
    return 0.0;
  }

  double rs = rs_of_density(n);
  if (rs < 1.0) {
    const double A = 0.0311, B = -0.048, C = 0.0020, D = -0.0116;

    return A * log(rs) + (B - A / 3.0) + (2.0 / 3.0) * C * rs * log(rs) +
           (2.0 * D - C) / 3.0 * rs;
  } else {
    const double gamma = -0.1423, beta1 = 1.0529, beta2 = 0.3334;
    double sq = sqrt(rs);
    double denom = 1.0 + beta1 * sq + beta2 * rs;
    double ec = gamma / denom;

    return ec * (1.0 + (7.0 / 6.0) * beta1 * sq + (4.0 / 3.0) * beta2 * rs) /
           denom;
  }
}

double lda_xc_energy_density(double n) {
  return lda_exchange_energy_density(n) +
         lda_correlation_energy_density_pz81(n);
}

double lda_xc_potential(double n) {
  return lda_exchange_potential(n) + lda_correlation_potential_pz81(n);
}

/* ---------------------------------------------------------------------
 * PW92 LDA correlation (\zeta=0), and PBE GGA exchange-correlation built on it.
 * Reference: Perdew & Wang, Phys. Rev. B 45, 13244 (1992); Perdew, Burke &
 * Ernzerhof (PBE), Phys. Rev. Lett. 77, 3865 (1996).
 * ------------------------------------------------------------------- */

/* PW92's G(rs) parametrization (unpolarized/zeta=0 branch) and its
 * rs-derivative, both needed (latter for LDA potential and as PBE correlation's
 * LDA-limit ingredient). Closed forms below (reference PBE subroutine CORPBE):
 *   Q1 = 2 * A * (b1 * \sqrt(rs) + b2 * rs + b3 * rs^(3/2) + b4 * rs^2),
 *   ec = -2 * A * (1+ a1 *rs) * \ln(1 + 1 / Q1),
 * and its derivative follows by direct product + chain rule differentiation of
 * that closed form.
 */
static void pw92_g_and_derivative(double rs, double *g, double *dg_drs) {
  const double A = 0.031091, a1 = 0.21370;
  const double b1 = 7.5957, b2 = 3.5876, b3 = 1.6382, b4 = 0.49294;

  double sq = sqrt(rs);
  double Q0 = -2.0 * A * (1.0 + a1 * rs);
  double Q1 = 2.0 * A * (b1 * sq + b2 * rs + b3 * rs * sq + b4 * rs * rs);
  double Q2 = log1p(1.0 / Q1);
  double Q3 = A * (b1 / sq + 2.0 * b2 + 3.0 * b3 * sq + 4.0 * b4 * rs);

  *g = Q0 * Q2;
  *dg_drs = -2.0 * A * a1 * Q2 - Q0 * Q3 / (Q1 * Q1 + Q1);
}

double pw92_correlation_energy_density(double n) {
  if (n <= 0.0) {
    return 0.0;
  }

  double rs = rs_of_density(n);
  double ec, decdrs;

  pw92_g_and_derivative(rs, &ec, &decdrs);

  return ec;
}

double pw92_correlation_potential(double n) {
  if (n <= 0.0) {
    return 0.0;
  }

  double rs = rs_of_density(n);
  double ec, decdrs;

  pw92_g_and_derivative(rs, &ec, &decdrs);

  // V_c = d(n * ec) / dn = ec + n * (dec / drs) * (drs / dn)
  //     = ec - (rs / 3) * decdrs,
  // since drs/dn = -rs / (3 * n).
  return ec - (rs / 3.0) * decdrs;
}

double pbe_exchange_energy_density(double n, double sigma) {
  if (n <= 0.0) {
    return 0.0;
  }

  const double kappa = 0.804, mu = 0.2195149727645171;

  double kF = cbrt(3.0 * M_PI * M_PI * n);
  double grad_n = sqrt(sigma);
  double s = grad_n / (2.0 * kF * n);
  double u = mu * s * s / kappa;
  double Fx = 1.0 + kappa - kappa / (1.0 + u);

  return lda_exchange_energy_density(n) * Fx;
}

void pbe_exchange_potential(double n, double sigma, double *vrho,
                            double *vsigma) {
  if (n <= 0.0) {
    *vrho = 0.0;
    *vsigma = 0.0;
    return;
  }

  const double kappa = 0.804, mu = 0.2195149727645171;

  double kF = cbrt(3.0 * M_PI * M_PI * n);
  double grad_n = sqrt(sigma);
  double s = grad_n / (2.0 * kF * n);
  double u = mu * s * s / kappa;
  double Fx = 1.0 + kappa - kappa / (1.0 + u);
  double dFx_ds = 2.0 * mu * s / ((1.0 + u) * (1.0 + u));

  //   s(n, \sigma) = \sqrt(\sigma) / (2 * kF(n) * n)
  // since kF * n scales as n^(4/3),
  //   d(s)/dn = -(4/3) * s / n
  // and d(s) / d\sigma = s / (2 * \sigma) (\sigma>0).
  double ds_dn = -(4.0 / 3.0) * s / n;
  double ds_dsigma = sigma > 0.0 ? s / (2.0 * sigma) : 0.0;

  double eps_x_lda = lda_exchange_energy_density(n);
  double vx_lda = lda_exchange_potential(n); // d(n*eps_x_lda)/dn

  *vrho = vx_lda * Fx + n * eps_x_lda * dFx_ds * ds_dn;
  *vsigma = n * eps_x_lda * dFx_ds * ds_dsigma;
}

double pbe_correlation_energy_density(double n, double sigma) {
  if (n <= 0.0) {
    return 0.0;
  }

  const double gamma = (1.0 - M_LN2) / (M_PI * M_PI);
  const double beta = 0.06672455060314922;

  double rs = rs_of_density(n);
  double ec;
  double decdrs;

  pw92_g_and_derivative(rs, &ec, &decdrs);

  double kF = cbrt(3.0 * M_PI * M_PI * n);
  double ks = sqrt(4.0 * kF / M_PI);
  double grad_n = sqrt(sigma);
  double t = grad_n / (2.0 * ks * n);
  double T2 = t * t;

  double y = exp(-ec / gamma);
  double B = (beta / gamma) / (y - 1.0);
  double Q4 = 1.0 + B * T2;
  double Q5 = 1.0 + B * T2 + B * B * T2 * T2;
  double X = (beta / gamma) * T2 * Q4 / Q5;
  double H = gamma * log1p(X);

  return ec + H;
}

void pbe_correlation_potential(double n, double sigma, double *vrho,
                               double *vsigma) {
  if (n <= 0.0) {
    *vrho = 0.0;
    *vsigma = 0.0;
    return;
  }

  const double gamma = (1.0 - M_LN2) / (M_PI * M_PI);
  const double beta = 0.06672455060314922;

  double rs = rs_of_density(n);
  double ec;
  double decdrs;

  pw92_g_and_derivative(rs, &ec, &decdrs);

  double kF = cbrt(3.0 * M_PI * M_PI * n);
  double ks = sqrt(4.0 * kF / M_PI);
  double grad_n = sqrt(sigma);
  double t = grad_n / (2.0 * ks * n);
  double T2 = t * t;

  //  B(ec) = (\beta / \gamma) / (\exp(-ec / \gamma) - 1)
  // and dB/dec by the quotient + chain rule
  double y = exp(-ec / gamma);
  double B = (beta / gamma) / (y - 1.0);
  double dB_dec = beta * y / (gamma * gamma * (y - 1.0) * (y - 1.0));

  double Q4 = 1.0 + B * T2;
  double Q5 = 1.0 + B * T2 + B * B * T2 * T2;
  double X = (beta / gamma) * T2 * Q4 / Q5;
  double H = gamma * log1p(X);

  // X = (beta/gamma) * f(T2,B)/Q5(T2,B), f = T2*Q4 = T2 + B*T2^2.
  double f = T2 * Q4;
  double df_dT2 = 1.0 + 2.0 * B * T2;
  double dQ5_dT2 = B + 2.0 * B * B * T2;
  double dX_dT2 = (beta / gamma) * (df_dT2 * Q5 - f * dQ5_dT2) / (Q5 * Q5);

  double df_dB = T2 * T2;
  double dQ5_dB = T2 + 2.0 * B * T2 * T2;
  double dX_dB = (beta / gamma) * (df_dB * Q5 - f * dQ5_dB) / (Q5 * Q5);

  double dH_dT2 = gamma / (1.0 + X) * dX_dT2;
  double dH_dB = gamma / (1.0 + X) * dX_dB;

  // Chain rule to n and \sigma:
  //  T2 = \sigma / (4 * ks(n)^2 * n^2) scales as n^(-7/3)
  // at fixed \sigma (ks^2 * n^2 ~ n^(7/3)),
  // so d(T2)/dn = -(7/3) * T2/n;
  // B depends on n only through ec(rs(n)), drs/dn = -rs/(3n)
  double drs_dn = -rs / (3.0 * n);
  double decdn = decdrs * drs_dn;

  double dT2_dn = -(7.0 / 3.0) * T2 / n;
  double dT2_dsigma = sigma > 0.0 ? T2 / sigma : 1.0 / (4.0 * ks * ks * n * n);

  double dH_dn = dH_dT2 * dT2_dn + dH_dB * dB_dec * decdn;
  double dH_dsigma = dH_dT2 * dT2_dsigma;

  *vrho = (ec + H) + n * (decdn + dH_dn);
  *vsigma = n * dH_dsigma;
}

double pbe_xc_energy_density(double n, double sigma) {
  return pbe_exchange_energy_density(n, sigma) +
         pbe_correlation_energy_density(n, sigma);
}

void pbe_xc_potential(double n, double sigma, double *vrho, double *vsigma) {
  double vrho_x, vsigma_x, vrho_c, vsigma_c;

  pbe_exchange_potential(n, sigma, &vrho_x, &vsigma_x);
  pbe_correlation_potential(n, sigma, &vrho_c, &vsigma_c);

  *vrho = vrho_x + vrho_c;
  *vsigma = vsigma_x + vsigma_c;
}

/* ---------------------------------------------------------------------
 * SCF machinery : structurally parallel to hartree_fock.c's, with dense Fock
 * exchange (K) matrix replaced by local V_xc(n(r)) diagonal potential.
 * ------------------------------------------------------------------- */

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

/* Kohn-Sham Hamiltonian matrix at the current density (encoded via u[],
 * n_orbitals doubly-occupied orbitals). Diagonal-only potential (kinetic
 * off-diagonals aside), unlike hartree_fock.c's dense exchange matrix.
 */
static cmatrix_t *build_ks_matrix(const double *r, int N, double dr, double Z,
                                  double **u, int n_orbitals) {
  double coeff = 0.5;
  double diag_factor = 2.0 * coeff / (dr * dr);
  double offdiag_factor = -coeff / (dr * dr);

  cmatrix_t *H = cmatrix_alloc(N, N);
  if (!H) {
    return NULL;
  }

  double *V_H = malloc((size_t)N * sizeof(double));
  double *n_dens = malloc((size_t)N * sizeof(double));
  if (!V_H || !n_dens) {
    free(V_H);
    free(n_dens);
    cmatrix_free(H);

    return NULL;
  }

  memset(V_H, 0, (size_t)N * sizeof(double));

  /* Classical Hartree potential: V_H(r) = sum_k occ_k * Y0_kk(r), occ=2
   * per doubly-occupied orbital -- identical formula/kernel to
   * hartree_fock.c's direct (J) term, reused verbatim via compute_Y0. */
  double *Y0 = malloc((size_t)N * sizeof(double));
  if (!Y0) {
    free(V_H);
    free(n_dens);
    cmatrix_free(H);

    return NULL;
  }

  for (int k = 0; k < n_orbitals; k++) {
    compute_Y0(r, N, dr, u[k], u[k], Y0);

    for (int i = 0; i < N; i++) {
      V_H[i] += 2.0 * Y0[i];
    }
  }

  free(Y0);

  // n(r) = \sum_k 2 * u_k(r)^2 / (4 * \pi * r^2): total electron density in
  // electrons/bohr^3, from all doubly-occupied orbitals.
  for (int i = 0; i < N; i++) {
    double r2 = r[i] * r[i];
    double u2_sum = 0.0;

    for (int k = 0; k < n_orbitals; k++) {
      u2_sum += u[k][i] * u[k][i];
    }

    n_dens[i] = (r2 > 0.0) ? 2.0 * u2_sum / (4.0 * M_PI * r2) : 0.0;
  }

  for (int i = 0; i < N; i++) {
    double Vnuc = (r[i] > 0.0) ? -Z / r[i] : 0.0;
    double Vxc = lda_xc_potential(n_dens[i]);

    CMAT(H, i, i) = c_real(diag_factor + Vnuc + V_H[i] + Vxc);

    if (i > 0) {
      CMAT(H, i, i - 1) = c_real(offdiag_factor);
    }

    if (i < N - 1) {
      CMAT(H, i, i + 1) = c_real(offdiag_factor);
    }
  }

  free(V_H);
  free(n_dens);

  /* Dirichlet boundaries: fully decouple endpoints so u(r_min)=u(r_max)=0
   * exactly
   * WARN: avoids -Z/r_min singularity acting as a spurious deep attractive well
   * pinned to a single grid point, which would otherwise dominate the low end
   * of spectrum.
   */
  for (int j = 0; j < N; j++) {
    CMAT(H, 0, j) = c_zero();
    CMAT(H, N - 1, j) = c_zero();
    CMAT(H, j, 0) = c_zero();
    CMAT(H, j, N - 1) = c_zero();
  }

  CMAT(H, 0, 0) = c_real(1e6 * diag_factor);
  CMAT(H, N - 1, N - 1) = c_real(1e6 * diag_factor);

  return H;
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

dft_result_t *dft_lda_atom_s_orbitals(const double *r, int N, double Z,
                                      int n_orbitals, double mix, double tol,
                                      int max_iter) {
  if (!r || N < 10 || Z <= 0.0 || n_orbitals < 1 || mix <= 0.0 || mix > 1.0 ||
      max_iter < 1) {
    return NULL;
  }

  double dr = r[1] - r[0];
  if (dr <= 0.0) {
    return NULL;
  }

  double **u = malloc((size_t)n_orbitals * sizeof(double *));
  if (!u) {
    return NULL;
  }

  for (int k = 0; k < n_orbitals; k++) {
    u[k] = malloc((size_t)N * sizeof(double));
  }

  // Initial guess: lowest n_orbitals bare-nucleus (l=0) eigenstates
  double Zc = Z;
  eigen_t *eig0 = central_potential_radial_solve((double *)r, N, 0, 1.0, 1.0,
                                                 V_coulomb, &Zc);
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

  double *orbital_energies = malloc((size_t)n_orbitals * sizeof(double));
  if (!orbital_energies) {
    free_orbital_arrays(u, n_orbitals);

    return NULL;
  }

  int converged = 0;
  int iter = 0;
  double eps_prev_sum = 0.0;

  for (; iter < max_iter; iter++) {
    cmatrix_t *H = build_ks_matrix(r, N, dr, Z, u, n_orbitals);
    if (!H) {
      free_orbital_arrays(u, n_orbitals);
      free(orbital_energies);

      return NULL;
    }

    eigen_t *eig = cmatrix_eigh_generic(H);
    cmatrix_free(H);
    if (!eig || eig->n < n_orbitals) {
      if (eig) {
        eigen_free(eig);
      }

      free_orbital_arrays(u, n_orbitals);
      free(orbital_energies);

      return NULL;
    }

    double eps_sum = 0.0;
    for (int k = 0; k < n_orbitals; k++) {
      eps_sum += eig->eigenvalues[k];
    }

    for (int k = 0; k < n_orbitals; k++) {
      double *u_new = malloc((size_t)N * sizeof(double));
      for (int i = 0; i < N; i++) {
        u_new[i] = CMAT(eig->eigenvectors, i, k).re;
      }

      normalize_u(u_new, N, dr);

      // keep a consistent sign (positive just passs inner boundary)
      if (u_new[1] < 0.0) {
        for (int i = 0; i < N; i++) {
          u_new[i] = -u_new[i];
        }
      }

      for (int i = 0; i < N; i++) {
        u[k][i] = mix * u_new[i] + (1.0 - mix) * u[k][i];
      }

      normalize_u(u[k], N, dr);
      free(u_new);

      orbital_energies[k] = eig->eigenvalues[k];
    }

    eigen_free(eig);

    if (iter > 0 && fabs(eps_sum - eps_prev_sum) < tol) {
      converged = 1;
      iter++;

      break;
    }

    eps_prev_sum = eps_sum;
  }

  // Final energetics at the converged density
  double *n_dens = malloc((size_t)N * sizeof(double));
  double *V_H = malloc((size_t)N * sizeof(double));
  double *Y0 = malloc((size_t)N * sizeof(double));
  memset(V_H, 0, (size_t)N * sizeof(double));
  for (int k = 0; k < n_orbitals; k++) {
    compute_Y0(r, N, dr, u[k], u[k], Y0);

    for (int i = 0; i < N; i++) {
      V_H[i] += 2.0 * Y0[i];
    }
  }

  free(Y0);

  for (int i = 0; i < N; i++) {
    double r2 = r[i] * r[i];
    double u2_sum = 0.0;

    for (int k = 0; k < n_orbitals; k++) {
      u2_sum += u[k][i] * u[k][i];
    }

    n_dens[i] = (r2 > 0.0) ? 2.0 * u2_sum / (4.0 * M_PI * r2) : 0.0;
  }

  /* E_H = \int u_tot(r)^2 * V_H(r) dr,
   * Where
   *  u_tot^2 = \sum_k u_k^2
   *  n(r) * 4 * \pi * r^2 = 2 * \sum_k u_k(r)^2, so this integral over dr is
   *  (1/2) * \int n(r) V_H(r) 4 * \pi * r^2 dr, Hartree energy.
   */
  double E_hartree = 0.0, E_xc = 0.0, Vxc_n_integral = 0.0;
  for (int i = 0; i < N; i++) {
    double u2_sum = 0.0;

    for (int k = 0; k < n_orbitals; k++) {
      u2_sum += u[k][i] * u[k][i];
    }

    E_hartree += u2_sum * V_H[i] * dr;
    E_xc += 2.0 * u2_sum * lda_xc_energy_density(n_dens[i]) * dr;
    Vxc_n_integral += 2.0 * u2_sum * lda_xc_potential(n_dens[i]) * dr;
  }

  free(V_H);
  free(n_dens);

  double eps_sum = 0.0;
  for (int k = 0; k < n_orbitals; k++) {
    eps_sum += orbital_energies[k];
  }
  double total_energy = 2.0 * eps_sum - E_hartree + E_xc - Vxc_n_integral;

  dft_result_t *res = malloc(sizeof(dft_result_t));
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
  res->E_hartree = E_hartree;
  res->E_xc = E_xc;
  res->iterations = iter;
  res->converged = converged;

  res->orbitals = malloc((size_t)n_orbitals * sizeof(cvector_t *));
  for (int k = 0; k < n_orbitals; k++) {
    res->orbitals[k] = cvector_alloc(N);

    for (int i = 0; i < N; i++) {
      res->orbitals[k]->data[i] = c_real(u[k][i]);
    }
  }

  free_orbital_arrays(u, n_orbitals);

  return res;
}

void dft_result_free(dft_result_t *res) {
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

  free(res);
}
