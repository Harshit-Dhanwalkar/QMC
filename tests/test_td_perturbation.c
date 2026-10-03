/*
 * Test: time-dependent perturbation theory (Dyson series).
 *
 * 1. Closed-form helpers: phase integral (incl. \Delta -> 0 limit and
 *    continuity across series/direct switch), constant, harmonic and
 *    Gaussian-pulse first-order amplitudes
 * 2. Numerical first order against those closed forms (constant, harmonic on
 *    and off resonance, Gaussian pulse incl. adiabatic-suppression law)
 * 3. Structure of series: c^(k) scales as exactly eps^k; with series
 *    truncated at order K error against exactly solvable two-level problem
 *    scales as eps^(K+1); Schrodinger-picture conversion agrees with exact
 *    solution
 * 4. A 4-level problem with a genuinely time-dependent, non-commuting V(t):
 *    error and norm deficit scaling against direct integration of full
 *    Schrodinger equation
 * 5. Fermi's golden rule: a dense band of final states under a constant
 *    perturbation accumulates probability at rate of fermi_golden_rate()
 *    from perturbation.h
 * 6. Invalid-input handling
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

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.10g expected=%.10g err=%.2e\n", label, got, expected,
         err);
  if (!(err <= tol)) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

static double cdiff(complex_t a, complex_t b) {
  return hypot(a.re - b.re, a.im - b.im);
}

static void put(cmatrix_t *m, int i, int j, double re, double im) {
  cmatrix_set(m, i, j, c_new(re, im));
}

static cvector_t *basis_state(int n, int idx) {
  cvector_t *v = cvector_alloc(n);
  if (v) {
    v->data[idx] = c_one();
  }

  return v;
}

// perturbations (callbacks)
typedef struct {
  double eps; // overall strength
  double Om;  // harmonic frequency (harmonic callback)
  double tau; // pulse width (Gaussian callback)
} drive_params_t;

static void constant_V(double t, void *p, cmatrix_t *V) {
  (void)t;

  const drive_params_t *q = p;
  put(V, 0, 1, q->eps, 0.0);
  put(V, 1, 0, q->eps, 0.0);
}

static void harmonic_V(double t, void *p, cmatrix_t *V) {
  const drive_params_t *q = p;
  double v = q->eps * cos(q->Om * t);

  put(V, 0, 1, v, 0.0);
  put(V, 1, 0, v, 0.0);
}

static void gaussian_V(double t, void *p, cmatrix_t *V) {
  const drive_params_t *q = p;
  double v = q->eps * exp(-t * t / (2.0 * q->tau * q->tau));

  put(V, 0, 1, v, 0.0);
  put(V, 1, 0, v, 0.0);
}

static void nan_V(double t, void *p, cmatrix_t *V) {
  (void)t;
  (void)p;

  put(V, 0, 0, NAN, 0.0);
}

// 4-level non-commuting, time-dependent perturbation
#define N4 4
typedef struct {
  double eps;
  cmatrix_t *A;
  cmatrix_t *B;
} four_level_t;

static void four_level_V(double t, void *p, cmatrix_t *V) {
  const four_level_t *q = p;
  double ca = q->eps * cos(0.7 * t);
  double cb = q->eps * exp(-(t - 1.5) * (t - 1.5));

  for (int i = 0; i < N4 * N4; i++) {
    V->data[i] = c_new(ca * q->A->data[i].re + cb * q->B->data[i].re,
                       ca * q->A->data[i].im + cb * q->B->data[i].im);
  }
}

static cmatrix_t *make_hermitian(double offset) {
  cmatrix_t *M = cmatrix_alloc(N4, N4);
  for (int i = 0; i < N4; i++) {
    for (int j = i; j < N4; j++) {
      double re = sin(offset + i + 2.0 * j);
      double im = (i < j) ? cos(offset + 2.0 + 3.0 * i + j) : 0.0;

      put(M, i, j, re, im);
      put(M, j, i, re, -im);
    }
  }

  return M;
}

// Direct RK4 integration of full Schrodinger equation (reference)
static void direct_evolve(const double *E, tdpt_perturbation_fn Vfn,
                          void *params, complex_t *psi, double t0, double t1,
                          int steps) {
  cmatrix_t *H0 = cmatrix_alloc(N4, N4);
  cmatrix_t *Hm = cmatrix_alloc(N4, N4);
  cmatrix_t *H1 = cmatrix_alloc(N4, N4);

  double h = (t1 - t0) / steps;
  complex_t k[4][N4];
  complex_t tmp[N4];
  for (int s = 0; s < steps; s++) {
    double t = t0 + s * h;

    memset(H0->data, 0, sizeof(complex_t) * N4 * N4);
    memset(Hm->data, 0, sizeof(complex_t) * N4 * N4);
    memset(H1->data, 0, sizeof(complex_t) * N4 * N4);

    Vfn(t, params, H0);
    Vfn(t + 0.5 * h, params, Hm);
    Vfn(t + h, params, H1);

    for (int i = 0; i < N4; i++) {
      H0->data[i * N4 + i].re += E[i];
      Hm->data[i * N4 + i].re += E[i];
      H1->data[i * N4 + i].re += E[i];
    }

    const cmatrix_t *Hs[4] = {H0, Hm, Hm, H1};
    const double coef[4] = {0.0, 0.5 * h, 0.5 * h, h};

    for (int stage = 0; stage < 4; stage++) {
      for (int i = 0; i < N4; i++) {
        complex_t v = psi[i];

        if (stage > 0) {
          v = c_add(v, c_scale(k[stage - 1][i], coef[stage]));
        }

        tmp[i] = v;
      }

      for (int i = 0; i < N4; i++) {
        complex_t acc = c_zero();

        for (int j = 0; j < N4; j++) {
          acc = c_add(acc, c_mul(Hs[stage]->data[i * N4 + j], tmp[j]));
        }

        k[stage][i] = c_new(acc.im, -acc.re); // -i * acc
      }
    }

    for (int i = 0; i < N4; i++) {
      complex_t inc = c_add(c_add(k[0][i], c_scale(k[1][i], 2.0)),
                            c_add(c_scale(k[2][i], 2.0), k[3][i]));
      psi[i] = c_add(psi[i], c_scale(inc, h / 6.0));
    }
  }

  cmatrix_free(H0);
  cmatrix_free(Hm);
  cmatrix_free(H1);
}

// Exact Schrodinger-picture solution of H = [[E0, V],[V, E1]] starting in |0>.
static void exact_two_level(double E0, double E1, double V, double T,
                            complex_t *a0, complex_t *a1) {
  double mean = 0.5 * (E0 + E1);
  double half = 0.5 * (E1 - E0);
  double W = sqrt(half * half + V * V);
  double c = cos(W * T);
  double s = sin(W * T);

  complex_t g = {cos(-mean * T), sin(-mean * T)};
  *a0 = c_mul(g, c_new(c, s * half / W)); // cos(WT) + i sin(WT) (Delta/2)/W
  *a1 = c_mul(g, c_new(0.0, -V * s / W));
}

// Tests
static void test_closed_forms(void) {
  printf("  === Closed-form helpers ===\n");

  complex_t z = tdpt_phase_integral(0.0, 3.0);
  check_close(z.re, 3.0, 1e-15, "phase integral at \\Delta = 0 (real)");
  check_close(z.im, 0.0, 1e-15, "phase integral at \\Delta = 0 (imag)");

  z = tdpt_phase_integral(1e-13, 3.0);
  check_close(z.re, 3.0, 1e-12, "tiny \\Delta, real part -> t");
  check_close(z.im, 4.5e-13, 1e-14, "tiny \\Delta, imag part -> \\Delta t^2/2");

  // continuity across the series/direct switch (x = \Delta t = 0.01)
  double t = 2.0;
  double below = 0.0099999999 / t;
  double above = 0.0100000001 / t;
  complex_t zb = tdpt_phase_integral(below, t);
  complex_t za = tdpt_phase_integral(above, t);
  double xb = below * t;
  complex_t direct = {sin(xb) / below, (1.0 - cos(xb)) / below};

  check_close(cdiff(zb, direct), 0.0, 1e-13, "series branch vs direct formula");
  check_close(cdiff(zb, za), 0.0, 1e-8, "continuous across the switch");

  // |\int|^2 = 4 \sin^2(x/2) / \Delta^2
  double d = 1.7;
  z = tdpt_phase_integral(d, 2.3);
  check_close(z.re * z.re + z.im * z.im,
              4.0 * pow(sin(0.5 * d * 2.3), 2) / (d * d), 1e-13,
              "|phase integral|^2");

  // constant perturbation: P = 4 V^2 \sin^2(w t / 2) / w^2
  complex_t c1 = tdpt_constant_first_order(c_new(0.03, 0.0), 1.3, 2.7);
  check_close(c1.re * c1.re + c1.im * c1.im,
              4.0 * 0.03 * 0.03 * pow(sin(0.5 * 1.3 * 2.7), 2) / (1.3 * 1.3),
              1e-14, "constant first-order probability (Rabi formula)");

  c1 = tdpt_constant_first_order(c_new(0.03, 0.0), 0.0, 2.7);
  check_close(c1.im, -0.03 * 2.7, 1e-14, "resonant constant: c = -i V t");

  // Gaussian adiabatic suppression law
  complex_t g0 = tdpt_gaussian_first_order(c_new(0.05, 0.0), 0.0, 0.9);
  complex_t g2 = tdpt_gaussian_first_order(c_new(0.05, 0.0), 2.0, 0.9);
  check_close(hypot(g2.re, g2.im) / hypot(g0.re, g0.im), exp(-0.5 * 4.0 * 0.81),
              1e-14, "Gaussian suppression exp(-w^2 tau^2 / 2)");

  complex_t bad = tdpt_gaussian_first_order(c_new(0.05, 0.0), 1.0, 0.0);
  check_true(isnan(bad.re) && isnan(bad.im), "Gaussian with tau = 0 -> NaN");

  bad = tdpt_phase_integral(NAN, 1.0);
  check_true(isnan(bad.re), "phase integral with NaN -> NaN");
}

static void test_first_order(void) {
  printf("  === Numerical first order vs closed forms ===\n");

  const double E[2] = {0.0, 1.3};
  cvector_t *psi0 = basis_state(2, 0);

  // constant coupling
  drive_params_t p = {0.05, 0.0, 0.0};
  tdpt_result_t *r = tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, 2.7, 400, 1);
  check_true(r != NULL, "constant: result returned");
  if (r) {
    complex_t ana = tdpt_constant_first_order(c_new(0.05, 0.0), 1.3, 2.7);
    check_close(cdiff(r->orders[1].data[1], ana), 0.0, 1e-9,
                "constant first order vs closed form");
    check_close(cdiff(r->orders[0].data[0], c_one()), 0.0, 1e-15,
                "order 0 is the initial state");
    check_close(cdiff(r->orders[1].data[0], c_zero()), 0.0, 1e-12,
                "no first-order change of the initial state (V_ii = 0)");

    tdpt_result_free(r);
  }

  // harmonic, off resonance and exactly resonant (\phi(0, t) = t limit)
  const double w_list[2] = {1.07, 1.0};
  for (int k = 0; k < 2; k++) {
    const double Eh[2] = {0.0, w_list[k]};
    drive_params_t ph = {0.04, 1.0, 0.0};

    r = tdpt_dyson(2, Eh, harmonic_V, &ph, psi0, 0.0, 6.0, 2000, 1);
    if (r) {
      complex_t ana =
          tdpt_harmonic_first_order(c_new(0.04, 0.0), w_list[k], 1.0, 6.0);
      char label[64];
      snprintf(label, sizeof label, "harmonic first order, \\omega_{fi} = %.2f",
               w_list[k]);

      check_close(cdiff(r->orders[1].data[1], ana), 0.0, 1e-10, label);

      tdpt_result_free(r);
    } else {
      check_true(0, "harmonic: result returned");
    }
  }

  // Gaussian pulse over the whole time axis
  const double Eg[2] = {0.0, 1.1};
  drive_params_t pg = {0.03, 0.0, 0.9};
  r = tdpt_dyson(2, Eg, gaussian_V, &pg, psi0, -9.0 * 0.9, 9.0 * 0.9, 4000, 1);
  if (r) {
    complex_t ana = tdpt_gaussian_first_order(c_new(0.03, 0.0), 1.1, 0.9);
    check_close(cdiff(r->orders[1].data[1], ana), 0.0, 1e-10,
                "Gaussian pulse first order vs closed form");

    tdpt_result_free(r);
  } else {
    check_true(0, "Gaussian: result returned");
  }

  cvector_free(psi0);
}

static void test_series_structure(void) {
  printf("  === Series structure (two-level, exactly solvable) ===\n");

  const double E[2] = {0.0, 1.3};
  const double T = 2.7;
  cvector_t *psi0 = basis_state(2, 0);

  // c^(k) is exactly homogeneous of degree k in V: halving eps scales it by
  // 2^-k
  drive_params_t pa = {0.04, 0.0, 0.0};
  drive_params_t pb = {0.02, 0.0, 0.0};
  tdpt_result_t *ra = tdpt_dyson(2, E, constant_V, &pa, psi0, 0.0, T, 400, 3);
  tdpt_result_t *rb = tdpt_dyson(2, E, constant_V, &pb, psi0, 0.0, T, 400, 3);
  if (ra && rb) {
    for (int k = 1; k <= 3; k++) {
      double na = cvector_norm(&ra->orders[k]);
      double nb = cvector_norm(&rb->orders[k]);
      char label[64];

      snprintf(label, sizeof label,
               "|c^(%d)| ratio for eps halving (expect %d)", k, 1 << k);

      check_close(na / nb, (double)(1 << k), 1e-8, label);
    }
  } else {
    check_true(0, "scaling runs returned results");
  }

  tdpt_result_free(ra);
  tdpt_result_free(rb);

  // truncation error ~ eps^(K+1) against the exact two-level solution
  for (int K = 1; K <= 3; K++) {
    double err[2];
    const double epss[2] = {0.08, 0.04};
    for (int j = 0; j < 2; j++) {
      drive_params_t p = {epss[j], 0.0, 0.0};
      tdpt_result_t *r = tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, T, 400, K);
      complex_t e0;
      complex_t e1;

      exact_two_level(E[0], E[1], epss[j], T, &e0, &e1);
      err[j] = INFINITY;

      if (r) {
        cvector_t *a = tdpt_schrodinger_amplitude(r, K);
        if (a) {
          err[j] = hypot(cdiff(a->data[0], e0), cdiff(a->data[1], e1));

          cvector_free(a);
        }
      }

      tdpt_result_free(r);
    }

    double ratio = err[0] / err[1];
    double expected = pow(2.0, K + 1);
    char label[96];
    snprintf(label, sizeof label,
             "K=%d truncation error ratio %.2f (expect %.0f)", K, ratio,
             expected);

    check_true(ratio > 0.8 * expected && ratio < 1.25 * expected, label);
  }

  // probabilities are picture independent and agree with the exact solution
  drive_params_t p = {0.03, 0.0, 0.0};
  tdpt_result_t *r = tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, T, 400, 3);
  if (r) {
    complex_t e0;
    complex_t e1;

    exact_two_level(E[0], E[1], 0.03, T, &e0, &e1);
    check_close(tdpt_probability(r, 3, 1), e1.re * e1.re + e1.im * e1.im, 1e-6,
                "P(1) through order 3 vs exact");
    check_close(tdpt_norm(r, 3), 1.0, 1e-6, "series norm through order 3");

    tdpt_result_free(r);
  }

  // max_order = 0 just returns the initial state
  r = tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, T, 10, 0);
  if (r) {
    check_close(cdiff(r->orders[0].data[0], c_one()), 0.0, 1e-15,
                "max_order 0 returns psi0");
    check_close(tdpt_norm(r, 0), 1.0, 1e-15, "max_order 0 norm");

    tdpt_result_free(r);
  } else {
    check_true(0, "max_order 0 returned a result");
  }

  cvector_free(psi0);
}

static void test_four_level(void) {
  printf("  === Four levels, time-dependent non-commuting V(t) ===");

  const double E[N4] = {0.0, 0.45, 1.1, 1.8};
  cmatrix_t *A = make_hermitian(0.3);
  cmatrix_t *B = make_hermitian(1.9);
  cvector_t *psi0 = cvector_alloc(N4);

  psi0->data[0] = c_new(0.5, 0.0);
  psi0->data[1] = c_new(0.0, 0.5);
  psi0->data[2] = c_new(0.5, 0.0);
  psi0->data[3] = c_new(-0.5, 0.0);

  const double T = 3.0;

  for (int K = 1; K <= 3; K++) {
    double err[2];
    double dev[2];
    const double epss[2] = {0.06, 0.03};
    for (int j = 0; j < 2; j++) {
      four_level_t fl = {epss[j], A, B};
      tdpt_result_t *r =
          tdpt_dyson(N4, E, four_level_V, &fl, psi0, 0.0, T, 600, K);

      // at t0 = 0 the interaction and Schrodinger pictures coincide
      complex_t ref[N4];
      memcpy(ref, psi0->data, sizeof ref);
      direct_evolve(E, four_level_V, &fl, ref, 0.0, T, 4000);

      err[j] = INFINITY;
      dev[j] = INFINITY;
      if (r) {
        cvector_t *a = tdpt_schrodinger_amplitude(r, K);
        if (a) {
          double s = 0.0;
          for (int i = 0; i < N4; i++) {
            double d = cdiff(a->data[i], ref[i]);
            s += d * d;
          }

          err[j] = sqrt(s);

          cvector_free(a);
        }

        dev[j] = fabs(tdpt_norm(r, K) - 1.0);
      }

      tdpt_result_free(r);
    }

    double ratio = err[0] / err[1];
    double expected = pow(2.0, K + 1);
    char label[96];

    snprintf(label, sizeof label,
             "K=%d error vs direct TDSE: ratio %.2f (expect %.0f)", K, ratio,
             expected);

    check_true(ratio > 0.8 * expected && ratio < 1.25 * expected, label);
    if (K <= 2) {
      double nratio = dev[0] / dev[1];
      snprintf(label, sizeof label,
               "K=%d norm deficit ratio %.2f (expect %.0f)", K, nratio,
               expected);

      check_true(nratio > 0.75 * expected && nratio < 1.3 * expected, label);
    }
  }

  cvector_free(psi0);
  cmatrix_free(A);
  cmatrix_free(B);
}

// Constant band-coupling perturbation: copies a fixed matrix.
static const cmatrix_t *g_band_matrix = NULL;

static void band_V(double t, void *p, cmatrix_t *V) {
  (void)t;
  (void)p;

  memcpy(V->data, g_band_matrix->data,
         (size_t)V->nrows * (size_t)V->ncols * sizeof(complex_t));
}

static void test_golden_rule(void) {
  printf("  === Fermi's golden rule from a dense band of final states ===\n");

  const int N = 200;
  const double dE = 0.02;
  const double Vc = 0.004;
  int n = N + 1;
  double *E = malloc((size_t)n * sizeof *E);
  cmatrix_t *Vm = cmatrix_alloc(n, n);
  cvector_t *psi0 = basis_state(n, 0);
  if (!E || !Vm || !psi0) {
    check_true(0, "allocation");

    free(E);
    cmatrix_free(Vm);
    cvector_free(psi0);
    return;
  }

  E[0] = 1.0;
  for (int j = 0; j < N; j++) {
    E[1 + j] = 1.0 + dE * (double)(j - N / 2);

    put(Vm, 1 + j, 0, Vc, 0.0);
    put(Vm, 0, 1 + j, Vc, 0.0);
  }

  double rate = fermi_golden_rate(Vm, 1, 0, 1.0 / dE);
  printf("  fermi_golden_rate = %.6e\n", rate);

  g_band_matrix = Vm;

  const double Ts[3] = {10.0, 20.0, 40.0};
  double ratio[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    tdpt_result_t *r = tdpt_dyson(n, E, band_V, NULL, psi0, 0.0, Ts[k],
                                  (int)(30.0 * Ts[k]), 1);
    if (r) {
      double P = 0.0;
      for (int f = 1; f < n; f++) {
        P += tdpt_probability(r, 1, f);
      }

      ratio[k] = P / (rate * Ts[k]);
      printf("  T=%4.0f  P/T = %.6e  ratio to golden rule = %.5f\n", Ts[k],
             P / Ts[k], ratio[k]);

      tdpt_result_free(r);
    } else {
      check_true(0, "band run returned a result");
    }
  }

  check_true(ratio[0] < ratio[1] && ratio[1] < ratio[2],
             "ratio approaches 1 monotonically as T grows");
  check_true(ratio[0] > 0.95 && ratio[0] < 1.0,
             "T = 10 within 5% (from below)");
  check_close(ratio[2], 1.0, 0.012, "T = 40 within 1.2% of the golden rule");

  free(E);
  cmatrix_free(Vm);
  cvector_free(psi0);
}

static void test_invalid(void) {
  printf("  === Invalid input ===\n");

  const double E[2] = {0.0, 1.0};
  const double Enan[2] = {0.0, NAN};
  drive_params_t p = {0.05, 0.0, 0.0};
  cvector_t *psi0 = basis_state(2, 0);
  cvector_t *wrong = cvector_alloc(3);
  cvector_t *nanpsi = basis_state(2, 0);

  nanpsi->data[1] = c_new(NAN, 0.0);

  check_true(tdpt_dyson(0, E, constant_V, &p, psi0, 0.0, 1.0, 10, 1) == NULL,
             "n = 0");
  check_true(tdpt_dyson(TDPT_MAX_DIM + 1, E, constant_V, &p, psi0, 0.0, 1.0, 10,
                        1) == NULL,
             "n above the limit");
  check_true(tdpt_dyson(2, NULL, constant_V, &p, psi0, 0.0, 1.0, 10, 1) == NULL,
             "NULL energies");
  check_true(tdpt_dyson(2, E, NULL, &p, psi0, 0.0, 1.0, 10, 1) == NULL,
             "NULL callback");
  check_true(tdpt_dyson(2, E, constant_V, &p, NULL, 0.0, 1.0, 10, 1) == NULL,
             "NULL psi0");
  check_true(tdpt_dyson(2, E, constant_V, &p, wrong, 0.0, 1.0, 10, 1) == NULL,
             "psi0 length mismatch");
  check_true(tdpt_dyson(2, E, constant_V, &p, nanpsi, 0.0, 1.0, 10, 1) == NULL,
             "non-finite psi0");
  check_true(tdpt_dyson(2, Enan, constant_V, &p, psi0, 0.0, 1.0, 10, 1) == NULL,
             "non-finite energy");
  check_true(tdpt_dyson(2, E, constant_V, &p, psi0, 1.0, 1.0, 10, 1) == NULL,
             "t1 == t0");
  check_true(tdpt_dyson(2, E, constant_V, &p, psi0, 1.0, 0.0, 10, 1) == NULL,
             "t1 < t0");
  check_true(tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, NAN, 10, 1) == NULL,
             "t1 = NaN");
  check_true(tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, 1.0, 3, 1) == NULL,
             "steps < 4");
  check_true(tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, 1.0, 10, -1) == NULL,
             "max_order < 0");
  check_true(tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, 1.0, 10,
                        TDPT_MAX_ORDER + 1) == NULL,
             "max_order above the limit");
  check_true(tdpt_dyson(2, E, nan_V, &p, psi0, 0.0, 1.0, 10, 1) == NULL,
             "callback producing NaN is rejected");

  tdpt_result_t *r = tdpt_dyson(2, E, constant_V, &p, psi0, 0.0, 1.0, 20, 2);
  check_true(r != NULL, "valid baseline run");
  check_true(tdpt_total_amplitude(NULL, 0) == NULL, "total_amplitude(NULL)");
  check_true(tdpt_total_amplitude(r, 3) == NULL, "order above max_order");
  check_true(tdpt_total_amplitude(r, -1) == NULL, "negative order");
  check_true(tdpt_schrodinger_amplitude(r, 3) == NULL,
             "schrodinger: bad order");
  check_true(isnan(tdpt_probability(r, 1, 2)),
             "probability: state out of range");
  check_true(isnan(tdpt_probability(r, 1, -1)), "probability: negative state");
  check_true(isnan(tdpt_probability(NULL, 0, 0)), "probability(NULL)");
  check_true(isnan(tdpt_norm(r, 5)), "norm: bad order");

  tdpt_result_free(r);
  tdpt_result_free(NULL);

  cvector_free(psi0);
  cvector_free(wrong);
  cvector_free(nanpsi);
}

int main(void) {
  printf(" > Time-dependent perturbation theory tests\n");

  test_closed_forms();
  test_first_order();
  test_series_structure();
  test_four_level();
  test_golden_rule();
  test_invalid();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll time-dependent perturbation checks passed\n");
  return 0;
}
