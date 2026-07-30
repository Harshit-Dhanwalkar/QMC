/*
Relativistic QM: Klein-Gordon and Dirac 1D solvers.
*/

#include "relativistic.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/linalg/linalg.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "potentials.h"
#include <math.h>
#include <stdlib.h>

eigen_t *klein_gordon_1d(double *x, int N, double *V, double m, double hbar,
                         double c) {
  if (!x || !V || N < 2) {
    return NULL;
  }

  // Discretize: (-\hbar^2 * c^2 d^2/dx^2 + m^2 * c^4 )\psi = ((E - V)^2)\psi
  //
  // The matrix below is the FREE (V-independent) operator -\hbar^2 c^2 d^2/dx^2
  // + m^2 c^4; it does not itself encode any spatial structure of V(x). Each
  // level n have expectation value :
  //   <V>_n = sum_i |c_i^(n)|^2 * V(x_i)
  // computed from level's (V-independent) eigenvector(i.e.1st-order
  // perturbative correction using free-particle-like modes as 0th-order basis).
  // Constant V(x) (every level reduces to <V>_n = V_avg identically, since
  // \sum_i |c_i|^2 = 1), and V(x) has real spatial structure
  //
  // NOTE: Know limitation : This is plain non-degenerate perturbation theory,
  // applied per level with no coupling between levels. For well-separated
  // levels it is real improvement over old uniform shift. But at
  // near-degeneracy, purely diagonal per-level correction can be worse than
  // uniform shift, because it's missing off-diagonal coupling that resolves
  // near-degenerate pair. On same test potential, ground state
  // (quasi-degenerate with first excited state) got worse under fix (~0.16 vs
  // ~0.05 Hartree-analog error vs. self-consistent reference) even though
  // levels 1,2,3,5 improved.
  // HACK:  If one need reliable ground-state accuracy on potential suspect has
  // close-lying levels, use klein_gordon_1d_self_consistent directly rather
  // than trusting this fast path's shift for that specific level.
  double dx = x[1] - x[0];
  double coeff = hbar * hbar * c * c / (dx * dx);
  double mc2 = m * c * c;

  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);

    return NULL;
  }

  for (int i = 0; i < N; i++) {
    diag[i] = 2.0 * coeff + mc2 * mc2;
    if (i < N - 1) {
      offdiag[i] = -coeff;
    }
  }

  // Need eigenvectors and eigenvalues to compute per-level <V>_n weighting
  eigen_t *eig = tridiag_eigh(diag, offdiag, N);
  free(diag);
  free(offdiag);
  if (!eig) {
    return NULL;
  }

  for (int n = 0; n < eig->n; n++) {
    double V_expect = 0.0;
    for (int i = 0; i < N; i++) {
      double c_i = CMAT(eig->eigenvectors, i, n).re;
      V_expect += c_i * c_i * V[i];
    }

    // E = <V>_n + \sqrt(\lambda_n) (positive-energy branch)
    eig->eigenvalues[n] = V_expect + sqrt(eig->eigenvalues[n]);
  }

  return eig;
}

klein_gordon_solution_t *
klein_gordon_1d_self_consistent(double *x, int N, double *V, double m,
                                double hbar, double c, double E_guess,
                                double tol, int max_iter) {
  if (!x || !V || N < 3 || max_iter < 1) {
    return NULL;
  }

  double dx = x[1] - x[0];
  double coeff = hbar * hbar * c * c / (dx * dx);
  double mc4 = m * m * c * c * c * c;
  int sign = (E_guess >= 0.0) ? 1 : -1;

  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);

    return NULL;
  }

  for (int i = 0; i < N - 1; i++) {
    offdiag[i] = -coeff;
  }

  double E = E_guess;
  double lambda_target = E * E;
  eigen_t *eig = NULL;
  int level = -1;
  int converged = 0;
  int iter;

  for (iter = 0; iter < max_iter; iter++) {
    // H(E)_diag[i] = 2 * coeff + m^2 * c^4 + V(x_i)^2 - 2 * E * V(x_i)
    for (int i = 0; i < N; i++) {
      diag[i] = 2.0 * coeff + mc4 + 2.0 * E * V[i] - V[i] * V[i];
    }

    if (eig) {
      eigen_free(eig);
    }
    eig = tridiag_eigh(diag, offdiag, N);
    if (!eig) {
      free(diag);
      free(offdiag);

      return NULL;
    }

    level = 0;
    double best = fabs(eig->eigenvalues[0] - lambda_target);
    for (int k = 1; k < eig->n; k++) {
      double diff = fabs(eig->eigenvalues[k] - lambda_target);
      if (diff < best) {
        best = diff;
        level = k;
      }
    }

    double lambda = eig->eigenvalues[level];
    if (lambda < 0.0) {
      // No real-energy solution at this iterate for tracked level
      break;
    }

    double E_new = sign * sqrt(lambda);
    if (fabs(E_new - E) < tol) {
      E = E_new;
      converged = 1;
      iter++;

      break;
    }

    E = E_new;
    lambda_target = E * E;
  }

  free(diag);
  free(offdiag);

  if (!eig || level < 0) {
    if (eig) {
      eigen_free(eig);
    }

    return NULL;
  }

  klein_gordon_solution_t *sol = malloc(sizeof *sol);
  if (!sol) {
    eigen_free(eig);

    return NULL;
  }

  cvector_t *psi = cvector_alloc(N);
  if (!psi) {
    free(sol);
    eigen_free(eig);

    return NULL;
  }

  for (int i = 0; i < N; i++) {
    psi->data[i] = CMAT(eig->eigenvectors, i, level);
  }
  eigen_free(eig);

  // dx-weighted normalization
  double norm = 0.0;
  for (int i = 0; i < N; i++) {
    norm += psi->data[i].re * psi->data[i].re;
  }
  norm *= dx;
  if (norm > 1e-300 && isfinite(norm)) {
    double inv = 1.0 / sqrt(norm);
    for (int i = 0; i < N; i++) {
      psi->data[i].re *= inv;
    }
  }

  sol->energy = E;
  sol->psi = psi;
  sol->converged = converged;
  sol->iterations = iter;

  return sol;
}

void klein_gordon_solution_free(klein_gordon_solution_t *sol) {
  if (!sol) {
    return;
  }
  cvector_free(sol->psi);
  free(sol);
}

eigen_t *dirac_1d(double *x, int N, double *V, double m, double hbar,
                  double c) {
  // Build 2N x 2N Hermitian matrix:
  // [ V(x) + mc^2,          -i * \hbar * c d/dx ]
  // [ -i * hbar * c * d/dx,  V(x) - mc^2        ]
  // via central differences.
  double dx = x[1] - x[0];
  int M = 2 * N;
  cmatrix_t *H = cmatrix_alloc(M, M);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

  double coeff = hbar * c / (2.0 * dx);
  for (int i = 0; i < N; i++) {
    int row1 = i, row2 = i + N;
    CMAT(H, row1, row1) = c_real(V[i] + m * c * c);
    CMAT(H, row2, row2) = c_real(V[i] - m * c * c);

    if (i > 0) {
      CMAT(H, row1, i - 1 + N) = c_imag(coeff);
      CMAT(H, row2, i - 1) = c_imag(coeff);
    }

    if (i < N - 1) {
      CMAT(H, row1, i + 1 + N) = c_imag(-coeff);
      CMAT(H, row2, i + 1) = c_imag(-coeff);
    }
  }

  eigen_t *eig = cmatrix_eigh_complex(H);
  cmatrix_free(H);

  return eig;
}

eigen_t *dirac_radial_solve(double *r, int N, int kappa, potential_fn V,
                            void *params, double m, double hbar, double c) {
  if (!r || N < 3 || !V || kappa == 0 || m <= 0.0 || hbar <= 0.0 || c <= 0.0) {
    return NULL;
  }

  double dr = r[1] - r[0];
  if (dr <= 0.0) {
    return NULL;
  }

  double mc2 = m * c * c;
  int M = 2 * N;
  cmatrix_t *H = cmatrix_alloc(M, M);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

  for (int i = 0; i < N; i++) {
    double Vi = V(r[i], params);
    CMAT(H, i, i) = c_real(Vi + mc2);
    CMAT(H, N + i, N + i) = c_real(Vi - mc2);

    double kc = hbar * c * kappa / r[i];
    CMAT(H, i, N + i) = c_add(CMAT(H, i, N + i), c_real(kc));
    CMAT(H, N + i, i) = c_add(CMAT(H, N + i, i), c_real(kc));

    if (i < N - 1) {
      double val = -hbar * c / (2.0 * dr);
      CMAT(H, i, N + i + 1) = c_add(CMAT(H, i, N + i + 1), c_real(val));
      CMAT(H, N + i + 1, i) = c_add(CMAT(H, N + i + 1, i), c_real(val));
    }
    if (i > 0) {
      double val = hbar * c / (2.0 * dr);
      CMAT(H, i, N + i - 1) = c_add(CMAT(H, i, N + i - 1), c_real(val));
      CMAT(H, N + i - 1, i) = c_add(CMAT(H, N + i - 1, i), c_real(val));
    }
  }

  eigen_t *eig = cmatrix_eigh_complex(H);
  cmatrix_free(H);

  return eig;
}

double dirac_hydrogen_energy_level(int n, int kappa, double Z, double hbar,
                                   double mass, double e_charge, double eps0,
                                   double c) {
  int abs_kappa = abs(kappa);
  if (abs_kappa < 1) {
    return NAN;
  }
  int n_r = n - abs_kappa; // radial quantum number
  if (kappa > 0) {
    if (n_r < 1) {
      return NAN; // l=\kappa branch: n_r=0 not allowed, need n > |\kappa|
    }
  } else {
    if (n_r < 0) {
      return NAN; // l=|\kappa|-1 branch: n_r=0 allowed, n = |\kappa| valid
                  // (ground state)
    }
  }

  double alpha = (e_charge * e_charge) / (4.0 * M_PI * eps0 * hbar * c);
  double Za = Z * alpha;
  double kappa2 = (double)(kappa * kappa);

  if (Za * Za >= kappa2) {
    return NAN; // beyond point-charge Dirac breakdown (Z * \alpha >= |\kappa|)
  }

  double denom = (double)(n - abs_kappa) + sqrt(kappa2 - Za * Za);
  double mc2 = mass * c * c;

  return mc2 / sqrt(1.0 + (Za * Za) / (denom * denom));
}
