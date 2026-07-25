/*
Numerov's method for Schrödinger equation

Reference: https://en.wikipedia.org/wiki/Numerov%27s_method
*/

#include "numerov.h"
#include "../../core/complex.h"
#include "../../core/linalg/tridiag_eigh.h"
#include "../../core/matrix.h"
#include "../../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 *   y_{n+1}*(1 - h^2 f_{n+1}/12) = 2*y_n*(1 + 5h^2 f_n/12)
 *                                  - y_{n-1}*(1 - h^2 f_{n-1}/12)
 */
void numerov_integrate(const numerov_params_t *p, double E, cvector_t *psi) {
  int N = p->n;
  double h2 = p->dx * p->dx;
  double *f = malloc(N * sizeof *f);
  if (!f) {
    return;
  }

  for (int i = 0; i < N; i++) {
    f[i] = (p->V[i] - E) / p->hbar_sq_2m;
  }

  psi->data[0].re = 0.0;
  psi->data[0].im = 0.0;
  psi->data[1].re = 1e-8;
  psi->data[1].im = 0.0;

  for (int i = 1; i < N - 1; i++) {
    double num = 2.0 * (1.0 + (5.0 / 12.0) * h2 * f[i]) * psi->data[i].re -
                 (1.0 - (1.0 / 12.0) * h2 * f[i - 1]) * psi->data[i - 1].re;
    double den = 1.0 - (1.0 / 12.0) * h2 * f[i + 1];
    double v = (fabs(den) > 1e-300) ? num / den : 0.0;
    if (!isfinite(v)) {
      v = 0.0;
    }

    psi->data[i + 1].re = v;
    psi->data[i + 1].im = 0.0;
    if (fabs(v) > 1e50) {
      for (int k = 0; k <= i + 1; k++) {
        psi->data[k].re *= 1e-50;
      }
    }
  }
  free(f);
}

// Find eigenstate using tridiagonal eigensolver
numerov_solution_t *numerov_shoot(numerov_params_t *params, double E_guess,
                                  double E_tol) {
  if (!params || !params->x || !params->V || params->n < 3) {
    return NULL;
  }

  int N = params->n;

  // Build tridiagonal Hamiltonian
  double coeff = params->hbar_sq_2m / (params->dx * params->dx);
  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);

    return NULL;
  }

  for (int i = 0; i < N; i++) {
    diag[i] = 2.0 * coeff + params->V[i];
  }

  for (int i = 0; i < N - 1; i++) {
    offdiag[i] = -coeff;
  }

  eigen_t *eig = tridiag_eigh(diag, offdiag, N);
  free(diag);
  free(offdiag);
  if (!eig)
    return NULL;

  int level = 0;
  double best_diff = fabs(eig->eigenvalues[0] - E_guess);
  for (int k = 1; k < eig->n; k++) {
    double diff = fabs(eig->eigenvalues[k] - E_guess);
    if (diff < best_diff) {
      best_diff = diff;
      level = k;
    }
  }

  // (void)E_tol;
  double E = eig->eigenvalues[level];
  cvector_t *psi = cvector_alloc(N);
  if (!psi) {
    eigen_free(eig);

    return NULL;
  }

  // Extract eigenvector column
  for (int i = 0; i < N; i++) {
    psi->data[i] = CMAT(eig->eigenvectors, i, level);
  }
  eigen_free(eig);

  // Normalize continuum convention
  double norm = 0.0;
  for (int i = 0; i < N; i++) {
    norm += psi->data[i].re * psi->data[i].re;
  }

  norm *= params->dx;
  if (norm > 1e-300 && isfinite(norm)) {
    double inv = 1.0 / sqrt(norm);
    for (int i = 0; i < N; i++) {
      psi->data[i].re *= inv;
      psi->data[i].im = 0.0;
    }
  }

  // NOTE: First significant lobe positive
  for (int i = 1; i < N; i++)
    if (fabs(psi->data[i].re) > 1e-10) {
      if (psi->data[i].re < 0.0) {
        for (int j = 0; j < N; j++) {
          psi->data[j].re = -psi->data[j].re;
        }
      }
      break;
    }

  numerov_solution_t *sol = malloc(sizeof *sol);
  if (!sol) {
    cvector_free(psi);

    return NULL;
  }

  sol->energy = E;
  sol->psi = psi;

  return sol;
}

void numerov_solution_free(numerov_solution_t *sol) {
  if (!sol) {
    return;
  }
  cvector_free(sol->psi);

  free(sol);
}

// Outer classical turning point
static int find_turning_point_index(const double *V, double E, int N) {
  int idx = N / 2;
  for (int i = N - 2; i > N / 2; i--) {
    if (V[i] < E && V[i + 1] >= E) {
      idx = i;

      break;
    }
  }

  int lo = N / 4, hi = 3 * N / 4;
  if (idx < lo) {
    idx = lo;
  }
  if (idx > hi) {
    idx = hi;
  }

  return idx;
}

// Log-derivative mismatch at the turning point for trial energy E.
static int log_deriv_mismatch(const numerov_params_t *params, double E,
                              double *mismatch, int *out_m_idx) {
  int N = params->n;
  int m_idx = find_turning_point_index(params->V, E, N);

  // Forward integration through turning point
  int n_left = m_idx + 2;
  numerov_params_t left_params = {.x = NULL,
                                  .V = params->V,
                                  .n = n_left,
                                  .dx = params->dx,
                                  .hbar_sq_2m = params->hbar_sq_2m};
  cvector_t *left = cvector_alloc(n_left);
  if (!left) {
    return 0;
  }
  numerov_integrate(&left_params, E, left);

  double y_left_m = left->data[m_idx].re;
  double y_left_p1 = left->data[m_idx + 1].re;
  double y_left_m1 = left->data[m_idx - 1].re;
  cvector_free(left);

  // Forward integration on reversed potential
  int n_right = N - m_idx + 1;
  double *V_rev = malloc(n_right * sizeof *V_rev);
  cvector_t *right_rev = cvector_alloc(n_right);
  if (!V_rev || !right_rev) {
    free(V_rev);
    cvector_free(right_rev);

    return 0;
  }

  for (int i = 0; i < n_right; i++) {
    V_rev[i] = params->V[N - 1 - i];
  }
  numerov_params_t right_params = {.x = NULL,
                                   .V = V_rev,
                                   .n = n_right,
                                   .dx = params->dx,
                                   .hbar_sq_2m = params->hbar_sq_2m};
  numerov_integrate(&right_params, E, right_rev);

  // Un-reverse indices: (m_idx-1, m_idx, m_idx+1)
  double right_m1 = right_rev->data[n_right - 1].re; // at m_idx-1
  double right_m = right_rev->data[n_right - 2].re;  // at m_idx
  double right_p1 = right_rev->data[n_right - 3].re; // at m_idx+1

  free(V_rev);
  cvector_free(right_rev);

  if (fabs(right_m) < 1e-300) {
    return 0;
  }

  // Rescale right segment
  double scale = y_left_m / right_m;
  double r_m1 = right_m1 * scale;
  double r_p1 = right_p1 * scale;

  double dL = (y_left_p1 - y_left_m1) / (2.0 * params->dx) / y_left_m;
  double dR = (r_p1 - r_m1) / (2.0 * params->dx) / y_left_m;

  *mismatch = dL - dR;
  *out_m_idx = m_idx;

  return 1;
}

numerov_solution_t *numerov_shoot_matching(numerov_params_t *params,
                                           double E_min, double E_max,
                                           int n_scan, double tol) {
  if (!params || !params->V || params->n < 8 || n_scan < 2 || E_max <= E_min) {
    return NULL;
  }

  double dE = (E_max - E_min) / (n_scan - 1);
  double prev_E = E_min, prev_val;
  int dummy_idx;
  if (!log_deriv_mismatch(params, prev_E, &prev_val, &dummy_idx)) {
    return NULL;
  }

  const double pole_thresh = 10.0;
  double lo = 0.0, hi = 0.0, m_lo = 0.0;
  int found = 0;

  for (int i = 1; i < n_scan && !found; i++) {
    double E = E_min + i * dE;
    double val;
    if (!log_deriv_mismatch(params, E, &val, &dummy_idx)) {
      prev_E = E;
      continue;
    }

    if (prev_val * val < 0.0 && fabs(prev_val) < pole_thresh &&
        fabs(val) < pole_thresh) {
      lo = prev_E;
      hi = E;
      m_lo = prev_val;
      found = 1;
    }

    prev_E = E;
    prev_val = val;
  }

  if (!found) {
    return NULL;
  }

  // Bisect within bracket to tolerance
  double E_final = 0.5 * (lo + hi);
  const int max_iter = 200;
  for (int iter = 0; iter < max_iter; iter++) {
    double mid = 0.5 * (lo + hi);
    double m_mid;
    int m_idx;

    if (!log_deriv_mismatch(params, mid, &m_mid, &m_idx)) {
      break;
    }

    E_final = mid;
    if (hi - lo < tol) {
      break;
    }

    if (m_lo * m_mid <= 0.0) {
      hi = mid;
    } else {
      lo = mid;
      m_lo = m_mid;
    }
  }

  // Build full wavefunction at converged energy
  int N = params->n;
  int m_idx = find_turning_point_index(params->V, E_final, N);
  int n_left = m_idx + 2;
  numerov_params_t left_params = {.x = NULL,
                                  .V = params->V,
                                  .n = n_left,
                                  .dx = params->dx,
                                  .hbar_sq_2m = params->hbar_sq_2m};
  cvector_t *left = cvector_alloc(n_left);
  if (!left) {
    return NULL;
  }
  numerov_integrate(&left_params, E_final, left);

  int n_right = N - m_idx + 1;
  double *V_rev = malloc(n_right * sizeof *V_rev);
  cvector_t *right_rev = cvector_alloc(n_right);
  if (!V_rev || !right_rev) {
    cvector_free(left);
    free(V_rev);
    cvector_free(right_rev);

    return NULL;
  }

  for (int i = 0; i < n_right; i++) {
    V_rev[i] = params->V[N - 1 - i];
  }
  numerov_params_t right_params = {.x = NULL,
                                   .V = V_rev,
                                   .n = n_right,
                                   .dx = params->dx,
                                   .hbar_sq_2m = params->hbar_sq_2m};
  numerov_integrate(&right_params, E_final, right_rev);

  double y_left_m = left->data[m_idx].re;
  double right_m = right_rev->data[n_right - 2].re;
  double scale = (fabs(right_m) > 1e-300) ? y_left_m / right_m : 0.0;

  cvector_t *psi = cvector_alloc(N);
  if (!psi) {
    cvector_free(left);
    free(V_rev);
    cvector_free(right_rev);

    return NULL;
  }

  for (int i = 0; i <= m_idx; i++) {
    psi->data[i].re = left->data[i].re;
    psi->data[i].im = 0.0;
  }

  // right_rev index (N-1-i) corresponds to original index i, for
  // i = m_idx + 1 .. N-1
  for (int i = m_idx + 1; i < N; i++) {
    int rev_idx = N - 1 - i; // position of original index i within right_rev
    psi->data[i].re = right_rev->data[rev_idx].re * scale;
    psi->data[i].im = 0.0;
  }
  cvector_free(left);
  free(V_rev);
  cvector_free(right_rev);

  // Continuum (dx-weighted) normalization
  double norm = 0.0;
  for (int i = 0; i < N; i++) {
    norm += psi->data[i].re * psi->data[i].re;
  }
  norm *= params->dx;

  if (norm > 1e-300 && isfinite(norm)) {
    double inv = 1.0 / sqrt(norm);
    for (int i = 0; i < N; i++) {
      psi->data[i].re *= inv;
    }
  }

  // 1st significant lobe positive
  for (int i = 1; i < N; i++)
    if (fabs(psi->data[i].re) > 1e-10) {
      if (psi->data[i].re < 0.0) {
        for (int j = 0; j < N; j++) {
          psi->data[j].re = -psi->data[j].re;
        }
      }
      break;
    }

  numerov_solution_t *sol = malloc(sizeof *sol);
  if (!sol) {
    cvector_free(psi);

    return NULL;
  }

  sol->energy = E_final;
  sol->psi = psi;

  return sol;
}
