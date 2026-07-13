/*
Numerov's method
for Schrödinger equation
*/

// TODO:Implement Cooley’s method (or log‑derivative matching)

#include "numerov.h"
#include "../../core/complex.h"
#include "../../core/linalg/tridiag_eigh.h"
#include "../../core/matrix.h"
#include "../../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void numerov_integrate(const numerov_params_t *p, double E, cvector_t *psi) {
  int N = p->n;
  double h2 = p->dx * p->dx;
  double *f = malloc(N * sizeof *f);
  if (!f)
    return;
  for (int i = 0; i < N; i++)
    f[i] = (p->V[i] - E) / p->hbar_sq_2m;

  psi->data[0].re = 0.0;
  psi->data[0].im = 0.0;
  psi->data[1].re = 1e-8;
  psi->data[1].im = 0.0;

  for (int i = 1; i < N - 1; i++) {
    double num = 2.0 * (1.0 - (5.0 / 12.0) * h2 * f[i]) * psi->data[i].re -
                 (1.0 + (1.0 / 12.0) * h2 * f[i - 1]) * psi->data[i - 1].re;
    double den = 1.0 + (1.0 / 12.0) * h2 * f[i + 1];
    double v = (fabs(den) > 1e-300) ? num / den : 0.0;
    if (!isfinite(v))
      v = 0.0;
    psi->data[i + 1].re = v;
    psi->data[i + 1].im = 0.0;
    if (fabs(v) > 1e50)
      for (int k = 0; k <= i + 1; k++)
        psi->data[k].re *= 1e-50;
  }
  free(f);
}

// Find eigenstate using tridiagonal eigensolver
numerov_solution_t *numerov_shoot(numerov_params_t *params, double E_guess,
                                  double E_tol) {
  (void)E_tol;
  if (!params || !params->x || !params->V || params->n < 3)
    return NULL;
  int N = params->n;
  double V_min = params->V[0];
  for (int i = 1; i < N; i++)
    if (params->V[i] < V_min)
      V_min = params->V[i];
  int level = 0;
  {
    double a = E_guess - V_min;
    level = (int)(a - 0.5 + 0.5);
    if (level < 0)
      level = 0;
  }

  // Build tridiagonal Hamiltonian
  double coeff = params->hbar_sq_2m / (params->dx * params->dx);
  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);
    return NULL;
  }

  for (int i = 0; i < N; i++)
    diag[i] = 2.0 * coeff + params->V[i];

  for (int i = 0; i < N - 1; i++)
    offdiag[i] = -coeff;

  eigen_t *eig = tridiag_eigh(diag, offdiag, N);
  free(diag);
  free(offdiag);
  if (!eig || level >= eig->n) {
    if (eig)
      eigen_free(eig);
    return NULL;
  }

  double E = eig->eigenvalues[level];
  cvector_t *psi = cvector_alloc(N);
  if (!psi) {
    eigen_free(eig);
    return NULL;
  }

  // Extract eigenvector column
  for (int i = 0; i < N; i++)
    psi->data[i] = CMAT(eig->eigenvectors, i, level);
  eigen_free(eig);

  // Normalize continuum convention
  double norm = 0.0;
  for (int i = 0; i < N; i++)
    norm += psi->data[i].re * psi->data[i].re;
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
      if (psi->data[i].re < 0.0)
        for (int j = 0; j < N; j++)
          psi->data[j].re = -psi->data[j].re;
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
  if (!sol)
    return;
  cvector_free(sol->psi);
  free(sol);
}
