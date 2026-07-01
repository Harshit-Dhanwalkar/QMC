/*
Numerov's method
for Schrödinger equation
*/

#include "numerov.h"
#include "../../core/complex.h"
#include "../../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Count zero crossings - equals the quantum number n for the n-th state
static int count_nodes(const cvector_t *psi, int N) {
  int nodes = 0;
  for (int i = 1; i < N; i++)
    if (psi->data[i - 1].re * psi->data[i].re < 0.0)
      nodes++;
  return nodes;
}

// Single forward Numerov integration at energy E
static double integrate_once(const numerov_params_t *p, double E,
                             cvector_t *psi) {
  int N = p->n;
  double h2 = p->dx * p->dx;

  double *f = malloc(N * sizeof *f);
  if (!f)
    return 0.0;

  for (int i = 0; i < N; i++)
    f[i] = (p->V[i] - E) / p->hbar_sq_2m;

  psi->data[0].re = 0.0;
  psi->data[0].im = 0.0;
  psi->data[1].re = 1e-6;
  psi->data[1].im = 0.0;

  for (int i = 1; i < N - 1; i++) {
    double num = 2.0 * (1.0 - (5.0 / 12.0) * h2 * f[i]) * psi->data[i].re -
                 (1.0 + (1.0 / 12.0) * h2 * f[i - 1]) * psi->data[i - 1].re;
    double denom = 1.0 + (1.0 / 12.0) * h2 * f[i + 1];
    psi->data[i + 1].re = (fabs(denom) > 1e-30) ? num / denom : 0.0;
    psi->data[i + 1].im = 0.0;
  }

  free(f);
  return psi->data[N - 1].re;
}

static int find_bracket(const numerov_params_t *p, int level, double *E_lo_out,
                        double *E_hi_out, cvector_t *psi) {
  double V_min = p->V[0];
  for (int i = 1; i < p->n; i++)
    if (p->V[i] < V_min)
      V_min = p->V[i];

  const double dE = 0.02;
  const double E_max = V_min + 500.0;

  double E_prev = V_min + 1e-8;
  double r_prev = integrate_once(p, E_prev, psi);

  for (double E = E_prev + dE; E < E_max; E += dE) {
    double r = integrate_once(p, E, psi);
    int nodes = count_nodes(psi, p->n);

    if (r_prev * r < 0.0 && nodes == level) {
      *E_lo_out = E_prev;
      *E_hi_out = E;
      return 1;
    }
    E_prev = E;
    r_prev = r;
  }
  return 0;
}

// Numerov integrator for given energy E
numerov_solution_t *numerov_shoot(numerov_params_t *params, double E_guess,
                                  double E_tol) {
  if (!params || !params->x || !params->V || params->n < 3)
    return NULL;

  int N = params->n;
  cvector_t *psi = cvector_alloc(N);
  if (!psi)
    return NULL;

  /* Determine target level from E_guess.
   * HO: E_n = n + 0.5, so level = round(E_guess - V_min - 0.5), clamped >= 0.
   */
  int level = 0;
  {
    double V_min = params->V[0];
    for (int i = 1; i < N; i++)
      if (params->V[i] < V_min)
        V_min = params->V[i];
    double above = E_guess - V_min;
    if (above > 1.0)
      level = (int)(above - 0.5 + 0.5);
    if (level < 0)
      level = 0;
  }

  double E_lo, E_hi;
  if (!find_bracket(params, level, &E_lo, &E_hi, psi)) {
    fprintf(stderr, "numerov_shoot: cannot bracket level %d (E_guess=%.4f)\n",
            level, E_guess);
    cvector_free(psi);
    return NULL;
  }

  double r_lo = integrate_once(params, E_lo, psi);
  double E_mid = E_lo;

  for (int iter = 0; iter < 300; iter++) {
    E_mid = 0.5 * (E_lo + E_hi);
    if (E_hi - E_lo < E_tol)
      break;

    double r_mid = integrate_once(params, E_mid, psi);
    double r_lo2 = integrate_once(params, E_lo, psi);

    if (r_lo2 * r_mid < 0.0)
      E_hi = E_mid;
    else
      E_lo = E_mid;
    r_lo = r_lo2;
  }
  (void)r_lo;

  integrate_once(params, E_mid, psi);

  // Normalize: sum|\psi_i|^2 * dx = 1 (continuum convention)
  double norm_sq = 0.0;
  for (int i = 0; i < N; i++)
    norm_sq += psi->data[i].re * psi->data[i].re;
  norm_sq *= params->dx;

  if (norm_sq > 1e-30) {
    double inv = 1.0 / sqrt(norm_sq);
    for (int i = 0; i < N; i++) {
      psi->data[i].re *= inv;
      psi->data[i].im = 0.0;
    }
  }

  // Sign convention: first significant lobe positive
  for (int i = 1; i < N; i++) {
    if (fabs(psi->data[i].re) > 1e-10) {
      if (psi->data[i].re < 0.0)
        for (int j = 0; j < N; j++)
          psi->data[j].re = -psi->data[j].re;
      break;
    }
  }

  numerov_solution_t *sol = malloc(sizeof *sol);
  if (!sol) {
    cvector_free(psi);
    return NULL;
  }
  sol->energy = E_mid;
  sol->psi = psi;
  return sol;
}

void numerov_solution_free(numerov_solution_t *sol) {
  if (!sol)
    return;
  cvector_free(sol->psi);
  free(sol);
}
