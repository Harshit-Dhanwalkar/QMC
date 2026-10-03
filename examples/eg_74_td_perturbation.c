/*
 * Time-dependent perturbation theory: Dyson series
 *
 * NOTE: eg_10_perturbation.c treats time-independent perturbations. Here
 * H(t) = H_0 + V(t) and system starts in an eigenstate of H_0. In interaction
 * picture amplitudes obey c' = -i V_I(t) c, and iterating this gives Dyson
 * series c = c^(0) + c^(1) + c^(2) + ..., where c^(k) is exactly k-th order in
 * V. All orders are integrated together, so any depth is as cheap as one extra
 * matrix-vector product per step
 *
 *  1. Where series works and where it fails: a two-level system with a constant
 *     coupling is exactly solvable, so each truncation order can be compared
 *     with truth. Near resonance V t stops being small and low orders lose
 *     oscillation
 *  2. Adiabatic suppression: a Gaussian pulse of width \tau transfers
 *     probability ~ \exp(-\omega^2 \tau^2); slow pulses leave system in its
 *     initial state
 *  3. Fermi's golden rule: coupling one state to a dense band of final states
 *     makes total probability grow linearly at rate 2 \pi |V|^2 \rho once t >>
 *     1/bandwidth
 */

#include "../physics/perturbation.h"
#include "../physics/td_perturbation.h"
#include "complex.h"
#include "matrix.h"
#include "vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void put(cmatrix_t *m, int i, int j, double re, double im) {
  cmatrix_set(m, i, j, c_new(re, im));
}

typedef struct {
  double eps;
  double tau;
} pulse_params_t;

static void constant_V(double t, void *p, cmatrix_t *V) {
  (void)t;

  const pulse_params_t *q = p;
  put(V, 0, 1, q->eps, 0.0);
  put(V, 1, 0, q->eps, 0.0);
}

static void gaussian_V(double t, void *p, cmatrix_t *V) {
  const pulse_params_t *q = p;
  double v = q->eps * exp(-t * t / (2.0 * q->tau * q->tau));

  put(V, 0, 1, v, 0.0);
  put(V, 1, 0, v, 0.0);
}

static const cmatrix_t *g_band = NULL;

static void band_V(double t, void *p, cmatrix_t *V) {
  (void)t;
  (void)p;

  memcpy(V->data, g_band->data,
         (size_t)V->nrows * (size_t)V->ncols * sizeof(complex_t));
}

static cvector_t *ground_state(int n) {
  cvector_t *v = cvector_alloc(n);
  if (v) {
    v->data[0] = c_one();
  }

  return v;
}

static void demo_series(void) {
  printf("  1. Series truncation vs exact two-level solution\n");
  printf("     E = {0, 0.2}, constant coupling V = 0.1 (V t ~ 1 at t ~ 10)\n");
  printf("     %-6s  %-10s  %-10s  %-10s  %-10s\n", "t", "exact P", "order 1",
         "order 2", "order 3");

  const double E[2] = {0.0, 0.2};
  const double V = 0.1;
  const double ts[5] = {1.0, 2.0, 5.0, 10.0, 20.0};
  cvector_t *psi0 = ground_state(2);
  pulse_params_t p = {V, 0.0};
  double half = 0.5 * (E[1] - E[0]);
  double W = sqrt(half * half + V * V);
  for (int k = 0; k < 5; k++) {
    tdpt_result_t *r = tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, ts[k],
                                  (int)(100 * ts[k]), 3);
    if (r) {
      double exact = (V / W) * (V / W) * pow(sin(W * ts[k]), 2);
      printf("     %-6.1f  %-10.6f  %-10.6f  %-10.6f  %-10.6f\n", ts[k], exact,
             tdpt_probability(r, 1, 1), tdpt_probability(r, 2, 1),
             tdpt_probability(r, 3, 1));

      tdpt_result_free(r);
    }
  }

  cvector_free(psi0);
  printf(
      "     (order 2 adds nothing to P(1): coupling is purely off-diagonal, so "
      "c^(2) only touches initial state. Orders agree while V t << 1, drift: "
      "near resonance at long times use Floquet or direct propagation)\n\n");
}

static void demo_adiabatic(void) {
  printf("  2. Adiabatic suppression by a Gaussian pulse (omega_fi = 1.1, V0 = "
         "0.03)\n");
  printf("     %-6s  %-14s  %-14s\n", "tau", "P numeric", "P closed form");
  const double E[2] = {0.0, 1.1};
  const double taus[4] = {0.5, 1.0, 2.0, 3.0};
  cvector_t *psi0 = ground_state(2);
  for (int k = 0; k < 4; k++) {
    pulse_params_t p = {0.03, taus[k]};
    tdpt_result_t *r = tdpt_dyson(2, E, gaussian_V, &p, psi0, -9.0 * taus[k],
                                  9.0 * taus[k], 4000, 1);
    if (r) {
      complex_t ana = tdpt_gaussian_first_order(c_new(0.03, 0.0), 1.1, taus[k]);
      printf("     %-6.1f  %-14.6e  %-14.6e\n", taus[k],
             tdpt_probability(r, 1, 1), ana.re * ana.re + ana.im * ana.im);

      tdpt_result_free(r);
    }
  }

  cvector_free(psi0);
  printf("     (probability falls like exp(-omega^2 tau^2): slow pulses are "
         "adiabatic)\n\n");
}

static void demo_golden_rule(void) {
  printf("  3. Fermi's golden rule from a band of 200 final states\n");
  const int N = 200;
  const double dE = 0.02;
  const double Vc = 0.004;
  int n = N + 1;
  double *E = malloc((size_t)n * sizeof *E);
  cmatrix_t *M = cmatrix_alloc(n, n);
  cvector_t *psi0 = ground_state(n);
  if (!E || !M || !psi0) {
    free(E);
    cmatrix_free(M);
    cvector_free(psi0);

    return;
  }

  E[0] = 1.0;
  for (int j = 0; j < N; j++) {
    E[1 + j] = 1.0 + dE * (double)(j - N / 2);

    put(M, 1 + j, 0, Vc, 0.0);
    put(M, 0, 1 + j, Vc, 0.0);
  }

  g_band = M;
  double rate = fermi_golden_rate(M, 1, 0, 1.0 / dE);
  printf("     golden-rule rate 2 pi V^2 rho = %.6e\n", rate);
  printf("     %-6s  %-14s  %-10s\n", "T", "P(T)/T", "ratio");

  const double Ts[3] = {10.0, 20.0, 40.0};
  for (int k = 0; k < 3; k++) {
    tdpt_result_t *r = tdpt_dyson(n, E, band_V, NULL, psi0, 0.0, Ts[k],
                                  (int)(30.0 * Ts[k]), 1);
    if (r) {
      double P = 0.0;

      for (int f = 1; f < n; f++) {
        P += tdpt_probability(r, 1, f);
      }

      printf("     %-6.0f  %-14.6e  %-10.4f\n", Ts[k], P / Ts[k],
             P / (rate * Ts[k]));

      tdpt_result_free(r);
    }
  }

  printf("     (ratio tends to 1 as T grows: rate emerges only once band is "
         "resolved, t >> 1/bandwidth)\n");

  free(E);
  cmatrix_free(M);
  cvector_free(psi0);
}

int main(void) {
  printf(" > Time-dependent perturbation theory (Dyson series)\n\n");

  demo_series();
  demo_adiabatic();
  demo_golden_rule();

  return 0;
}
