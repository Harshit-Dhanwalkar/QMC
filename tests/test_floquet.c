/*
 * Test: Floquet theory (time-domain and Sambe-space solvers)
 *
 * 1. J_0 helper against reference values
 * 2. Circular-drive Rabi problem (exactly solvable without RWA): both solvers
 *    reproduce quasi-energies \omega/2 +- \Omega_R/2 (mod omega); resonant
 *    quasi-energy gap equals \Omega; unitarity/eigen residuals are tiny
 * 3. RK4 convergence: halving step shrinks time-domain error ~16x
 * 4. Coherent destruction of tunneling: splitting = Delta |J_0(A/omega)| for
 *    omega >> Delta, vanishing at J_0 zero; time-domain and Sambe agree to
 *    ~1e-8 with each other
 * 5. Bloch-Siegert shift: resonance of linearly polarized drive sits at
 *    \omega_0 + \Omega^2/(4 \omega_0) to O(\Omega^4)
 * 6. Spin-1 with a bare spectrum wider than omega (folding across zone):
 *    solvers agree; stroboscopic evolution matches direct integration;
 *    modes are orthonormal; Sambe truncation weight falls with K
 * 7. Static and degenerate limits: folding of a static spectrum, H = 0
 * 8. Invalid-input handling, including non-Hermitian harmonics
 */

#include "../core/matrix.h"
#include "../physics/floquet.h"
#include "../physics/variational.h"
#include "complex.h"
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

static void put(cmatrix_t *m, int i, int j, double re, double im) {
  cmatrix_set(m, i, j, c_new(re, im));
}

// distance between two quasi-energies on circle of circumference \omega
static double circ_dist(double a, double b, double omega) {
  double d = fmod(fabs(a - b), omega);

  return fmin(d, omega - d);
}

// max |V^dagger V - I| for an n x n column-orthonormality check
static double orthonormality_error(const cmatrix_t *V) {
  int n = V->nrows;
  double err = 0.0;
  for (int a = 0; a < n; a++) {
    for (int b = 0; b < n; b++) {
      complex_t s = c_zero();

      for (int i = 0; i < n; i++) {
        s = c_add(s, c_mul(c_conj(V->data[i * n + a]), V->data[i * n + b]));
      }
      if (a == b) {
        s.re -= 1.0;
      }

      err = fmax(err, hypot(s.re, s.im));
    }
  }

  return err;
}

/* Hamiltonians (callbacks) */
typedef struct {
  double w0, Om, om;
} circ_params_t;

static void circular_H(double t, void *p, cmatrix_t *H) {
  const circ_params_t *q = p;

  put(H, 0, 0, 0.5 * q->w0, 0.0);
  put(H, 1, 1, -0.5 * q->w0, 0.0);
  put(H, 0, 1, 0.5 * q->Om * cos(q->om * t), -0.5 * q->Om * sin(q->om * t));
  put(H, 1, 0, 0.5 * q->Om * cos(q->om * t), 0.5 * q->Om * sin(q->om * t));
}

// Sambe harmonics (m = -1, 0, +1) of circular drive; caller frees
static void circular_harmonics(const circ_params_t *q, cmatrix_t *h[3]) {
  for (int i = 0; i < 3; i++) {
    h[i] = cmatrix_alloc(2, 2);
  }

  put(h[1], 0, 0, 0.5 * q->w0, 0.0);
  put(h[1], 1, 1, -0.5 * q->w0, 0.0);
  put(h[2], 1, 0, 0.5 * q->Om,
      0.0); // \exp^{+i w t} coefficient: (Om/2) \sigma_-
  put(h[0], 0, 1, 0.5 * q->Om,
      0.0); // \exp^{-i w t} coefficient: (Om/2) \sigma_+
}

static void free_harmonics(cmatrix_t **h, int count) {
  for (int i = 0; i < count; i++) {
    cmatrix_free(h[i]);
  }
}

typedef struct {
  double Delta, A, om;
} cdt_params_t;

static void cdt_H(double t, void *p, cmatrix_t *H) {
  const cdt_params_t *q = p;
  double z = 0.5 * q->A * cos(q->om * t);

  put(H, 0, 0, z, 0.0);
  put(H, 1, 1, -z, 0.0);
  put(H, 0, 1, 0.5 * q->Delta, 0.0);
  put(H, 1, 0, 0.5 * q->Delta, 0.0);
}

static void cdt_harmonics(const cdt_params_t *q, cmatrix_t *h[3]) {
  for (int i = 0; i < 3; i++) {
    h[i] = cmatrix_alloc(2, 2);
  }

  put(h[1], 0, 1, 0.5 * q->Delta, 0.0);
  put(h[1], 1, 0, 0.5 * q->Delta, 0.0);
  for (int i = 0; i < 3; i += 2) {
    put(h[i], 0, 0, 0.25 * q->A, 0.0);
    put(h[i], 1, 1, -0.25 * q->A, 0.0);
  }
}

typedef struct {
  double b, g, om;
} spin1_params_t;

static void spin1_H(double t, void *p, cmatrix_t *H) {
  const spin1_params_t *q = p;
  double c = q->g * cos(q->om * t) / sqrt(2.0);

  put(H, 0, 0, q->b, 0.0);
  put(H, 2, 2, -q->b, 0.0);
  put(H, 0, 1, c, 0.0);
  put(H, 1, 0, c, 0.0);
  put(H, 1, 2, c, 0.0);
  put(H, 2, 1, c, 0.0);
}

static void spin1_harmonics(const spin1_params_t *q, cmatrix_t *h[3]) {
  for (int i = 0; i < 3; i++) {
    h[i] = cmatrix_alloc(3, 3);
  }

  put(h[1], 0, 0, q->b, 0.0);
  put(h[1], 2, 2, -q->b, 0.0);

  double c = 0.5 * q->g / sqrt(2.0);
  for (int i = 0; i < 3; i += 2) {
    put(h[i], 0, 1, c, 0.0);
    put(h[i], 1, 0, c, 0.0);
    put(h[i], 1, 2, c, 0.0);
    put(h[i], 2, 1, c, 0.0);
  }
}

typedef struct {
  double eps0, eps1;
} static_params_t;

static void static_H(double t, void *p, cmatrix_t *H) {
  (void)t;
  const static_params_t *q = p;

  put(H, 0, 0, q->eps0, 0.0);
  put(H, 1, 1, q->eps1, 0.0);
}

static void zero_H(double t, void *p, cmatrix_t *H) {
  (void)t;
  (void)p;
  (void)H; // leaves H identically zero
}

/* Tests */
static void test_bessel(void) {
  printf("  === J0 helper ===\n");

  check_close(floquet_bessel_j0(0.0), 1.0, 1e-15, "J0(0)");
  check_close(floquet_bessel_j0(1.0), 0.7651976865579666, 1e-14, "J0(1)");
  check_close(floquet_bessel_j0(2.404825557695773), 0.0, 1e-14,
              "J0(first zero)");
  check_close(floquet_bessel_j0(10.0), -0.2459357644513483, 1e-14, "J0(10)");
  check_close(floquet_bessel_j0(-10.0), -0.2459357644513483, 1e-14, "J0 even");
  check_close(floquet_bessel_j0(50.0), 0.05581232766925181, 1e-13, "J0(50)");
  check_true(isnan(floquet_bessel_j0(NAN)), "J0(NaN) = NaN");
  check_true(isnan(floquet_bessel_j0(1.0e6)), "J0 out of range = NaN");
}

static void test_circular_rabi(void) {
  printf("  === Circular-drive Rabi problem (exact, no RWA) ===\n");

  circ_params_t p = {1.0, 0.3, 0.85};
  floquet_result_t *rt = floquet_solve_time(2, circular_H, &p, p.om, 2000);
  check_true(rt != NULL, "time solver returned a result");
  if (!rt) {
    return;
  }

  double lo = floquet_circular_rabi_quasienergy(p.w0, p.Om, p.om, -1);
  double hi = floquet_circular_rabi_quasienergy(p.w0, p.Om, p.om, +1);
  if (lo > hi) {
    double t = lo;
    lo = hi;
    hi = t;
  }

  check_close(rt->quasienergies[0], lo, 1e-9, "time: lower quasi-energy");
  check_close(rt->quasienergies[1], hi, 1e-9, "time: upper quasi-energy");
  check_true(rt->unitarity_error < 1e-10, "time: U(T) unitary");
  check_true(rt->eigen_residual < 1e-10, "time: eigen residual small");
  check_true(rt->quasienergies[0] >= -0.5 * p.om &&
                 rt->quasienergies[1] < 0.5 * p.om,
             "quasi-energies lie in Brillouin zone");
  check_close(orthonormality_error(rt->modes), 0.0, 1e-10,
              "time: modes orthonormal");

  cmatrix_t *h[3];
  circular_harmonics(&p, h);
  floquet_result_t *rs = floquet_solve_sambe(2, h, 1, p.om, 8);
  check_true(rs != NULL, "Sambe solver returned a result");
  if (rs) {
    check_close(rs->quasienergies[0], lo, 1e-12, "Sambe: lower quasi-energy");
    check_close(rs->quasienergies[1], hi, 1e-12, "Sambe: upper quasi-energy");
    check_true(rs->truncation_weight < 1e-12, "Sambe: negligible edge weight");

    floquet_result_free(rs);
  }

  free_harmonics(h, 3);
  floquet_result_free(rt);

  // Exact resonance: avoided-crossing gap is Rabi coupling \Omega
  circ_params_t pr = {1.0, 0.3, 1.0};
  floquet_result_t *rr = floquet_solve_time(2, circular_H, &pr, pr.om, 2000);
  if (rr) {
    check_close(circ_dist(rr->quasienergies[0], rr->quasienergies[1], pr.om),
                pr.Om, 1e-9, "resonant gap = Omega");

    floquet_result_free(rr);
  }
}

static void test_rk4_order(void) {
  printf("  === RK4 convergence order ===\n");

  circ_params_t p = {1.0, 0.3, 0.85};
  double exact_hi = floquet_circular_rabi_quasienergy(p.w0, p.Om, p.om, +1);
  double exact_lo = floquet_circular_rabi_quasienergy(p.w0, p.Om, p.om, -1);
  if (exact_lo > exact_hi) {
    exact_hi = exact_lo;
  }

  double err[2];
  const int steps[2] = {32, 64};
  for (int k = 0; k < 2; k++) {
    floquet_result_t *r = floquet_solve_time(2, circular_H, &p, p.om, steps[k]);
    err[k] = r ? fabs(r->quasienergies[1] - exact_hi) : INFINITY;

    printf("  steps=%d: quasi-energy error %.3e\n", steps[k], err[k]);

    floquet_result_free(r);
  }

  double ratio = err[0] / err[1];
  printf("  error ratio (expect ~16): %.2f\n", ratio);

  check_true(ratio > 10.0 && ratio < 24.0, "fourth-order convergence");
}

static void test_cdt(void) {
  printf("  === Coherent destruction of tunneling ===\n");

  double Delta = 0.1;
  double om = 4.0;
  const double xs[4] = {0.5, 1.5, 2.404825557695773, 3.5};

  for (int k = 0; k < 4; k++) {
    cdt_params_t p = {Delta, xs[k] * om, om};
    floquet_result_t *rt = floquet_solve_time(2, cdt_H, &p, om, 4000);

    cmatrix_t *h[3];
    cdt_harmonics(&p, h);
    floquet_result_t *rs = floquet_solve_sambe(2, h, 1, om, 30);
    if (rt && rs) {
      char label[96];
      double st = circ_dist(rt->quasienergies[0], rt->quasienergies[1], om);
      double ss = circ_dist(rs->quasienergies[0], rs->quasienergies[1], om);
      double pred = floquet_cdt_splitting(Delta, p.A, om);

      snprintf(label, sizeof label, "A/\\omega=%.4f time vs Sambe", xs[k]);
      check_close(st, ss, 1e-8, label);

      snprintf(label, sizeof label, "A/\\omega=%.4f Sambe vs \\Delta|J0|",
               xs[k]);
      check_close(ss, pred, 3e-5, label);
    } else {
      check_true(0, "CDT solvers returned results");
    }

    floquet_result_free(rt);
    floquet_result_free(rs);
    free_harmonics(h, 3);
  }

  check_close(floquet_cdt_splitting(0.1, 2.404825557695773 * 4.0, 4.0), 0.0,
              1e-14, "predicted splitting vanishes at J0 zero");
  check_true(isnan(floquet_cdt_splitting(0.1, 1.0, 0.0)),
             "splitting with \\omega = 0 is NaN");
}

typedef struct {
  double w0, Om;
  int K;
} bs_params_t;

static double bs_gap(double om, void *pp) {
  const bs_params_t *q = pp;
  cmatrix_t *h[3];
  for (int i = 0; i < 3; i++) {
    h[i] = cmatrix_alloc(2, 2);
  }

  put(h[1], 0, 0, 0.5 * q->w0, 0.0);
  put(h[1], 1, 1, -0.5 * q->w0, 0.0);
  for (int i = 0; i < 3; i += 2) {
    put(h[i], 0, 1, 0.5 * q->Om, 0.0);
    put(h[i], 1, 0, 0.5 * q->Om, 0.0);
  }

  floquet_result_t *r = floquet_solve_sambe(2, h, 1, om, q->K);
  double gap = INFINITY;
  if (r) {
    gap = circ_dist(r->quasienergies[0], r->quasienergies[1], om);
  }

  floquet_result_free(r);
  free_harmonics(h, 3);

  return gap;
}

static void test_bloch_siegert(void) {
  printf(" === Bloch-Siegert shift ===\n");

  const double Om_list[2] = {0.05, 0.1};
  for (int k = 0; k < 2; k++) {
    bs_params_t q = {1.0, Om_list[k], 6};
    double res = golden_section_minimize(0.99, 1.02, bs_gap, &q, 1e-10);
    double shift = res - q.w0;
    double pred = Om_list[k] * Om_list[k] / (4.0 * q.w0);
    char label[96];
    snprintf(label, sizeof label,
             "\\Omega=%.2f resonance shift vs \\Omega^2/(4w0)", Om_list[k]);

    // Next-order correction is O(\Omega^4/w0^3): ~1% at \Omega = 0.1
    check_true(fabs(shift - pred) < 0.02 * pred, label);
    printf("    shift=%.6e predicted=%.6e\n", shift, pred);
    snprintf(label, sizeof label, "\\Omega=%.2f minimum gap ~ \\Omega",
             Om_list[k]);
    check_close(bs_gap(res, &q), Om_list[k], 1e-3 * Om_list[k] + 1e-6, label);
  }
}

// Reference ODE integration of psi(t) with RK4 (independent of solvers).
static void direct_evolve(int n, floquet_hamiltonian_fn Hfn, void *params,
                          cvector_t *psi, double t_total, int steps) {
  cmatrix_t *H0 = cmatrix_alloc(n, n);
  cmatrix_t *Hm = cmatrix_alloc(n, n);
  cmatrix_t *H1 = cmatrix_alloc(n, n);

  double h = t_total / steps;
  cvector_t *k1 = cvector_alloc(n);
  cvector_t *k2 = cvector_alloc(n);
  cvector_t *k3 = cvector_alloc(n);
  cvector_t *k4 = cvector_alloc(n);
  cvector_t *tmp = cvector_alloc(n);
  for (int s = 0; s < steps; s++) {
    double t = s * h;
    memset(H0->data, 0, (size_t)n * n * sizeof(complex_t));
    memset(Hm->data, 0, (size_t)n * n * sizeof(complex_t));
    memset(H1->data, 0, (size_t)n * n * sizeof(complex_t));

    Hfn(t, params, H0);
    Hfn(t + 0.5 * h, params, Hm);
    Hfn(t + h, params, H1);

    const cmatrix_t *const Hs[4] = {H0, Hm, Hm, H1};
    cvector_t *ks[4] = {k1, k2, k3, k4};
    const double coef[4] = {0.0, 0.5 * h, 0.5 * h, h};
    for (int stage = 0; stage < 4; stage++) {
      for (int i = 0; i < n; i++) {
        complex_t v = psi->data[i];
        if (stage > 0) {
          v = c_add(v, c_scale(ks[stage - 1]->data[i], coef[stage]));
        }

        tmp->data[i] = v;
      }

      for (int i = 0; i < n; i++) {
        complex_t acc = c_zero();

        for (int j = 0; j < n; j++) {
          acc = c_add(acc, c_mul(Hs[stage]->data[i * n + j], tmp->data[j]));
        }

        ks[stage]->data[i] = c_new(acc.im, -acc.re); // -i * acc
      }
    }

    for (int i = 0; i < n; i++) {
      complex_t inc = c_add(c_add(k1->data[i], c_scale(k2->data[i], 2.0)),
                            c_add(c_scale(k3->data[i], 2.0), k4->data[i]));

      psi->data[i] = c_add(psi->data[i], c_scale(inc, h / 6.0));
    }
  }

  cmatrix_free(H0);
  cmatrix_free(Hm);
  cmatrix_free(H1);
  cvector_free(k1);
  cvector_free(k2);
  cvector_free(k3);
  cvector_free(k4);
  cvector_free(tmp);
}

static void test_spin1(void) {
  printf("  === Spin-1: bare spectrum wider than \\omega ===\n");

  spin1_params_t p = {2.3, 0.4, 1.1};
  floquet_result_t *rt = floquet_solve_time(3, spin1_H, &p, p.om, 4000);
  cmatrix_t *h[3];
  spin1_harmonics(&p, h);
  floquet_result_t *rs = floquet_solve_sambe(3, h, 1, p.om, 14);

  check_true(rt && rs, "both solvers returned results");
  if (!rt || !rs) {
    floquet_result_free(rt);
    floquet_result_free(rs);
    free_harmonics(h, 3);

    return;
  }

  for (int a = 0; a < 3; a++) {
    char label[64];
    snprintf(label, sizeof label, "quasi-energy %d: time vs Sambe", a);

    check_close(rt->quasienergies[a], rs->quasienergies[a], 1e-8, label);
  }

  check_true(rs->truncation_weight < 1e-12, "Sambe edge weight negligible");
  check_close(orthonormality_error(rt->modes), 0.0, 1e-10,
              "time modes orthonormal");
  check_close(orthonormality_error(rs->modes), 0.0, 1e-8,
              "Sambe modes orthonormal");

  // Spin-1 has symmetric spectrum {-\exp, 0, +\exp}
  check_close(rt->quasienergies[1], 0.0, 1e-8, "middle quasi-energy is 0");
  check_close(rt->quasienergies[0] + rt->quasienergies[2], 0.0, 1e-8,
              "spectrum symmetric");

  // Stroboscopic evolution against direct integration
  cvector_t *psi0 = cvector_alloc(3);
  psi0->data[0] = c_new(0.6, 0.0);
  psi0->data[1] = c_new(0.0, 0.48);
  psi0->data[2] = c_new(0.64, 0.0);

  cvector_t *direct = cvector_alloc(3);
  memcpy(direct->data, psi0->data, 3 * sizeof(complex_t));
  long nper = 7;
  direct_evolve(3, spin1_H, &p, direct, (double)nper * 2.0 * M_PI / p.om,
                (int)nper * 4000);
  cvector_t *strobe_t = floquet_stroboscopic_state(rt, psi0, nper);
  cvector_t *strobe_s = floquet_stroboscopic_state(rs, psi0, nper);
  if (strobe_t && strobe_s) {
    double dt = 0.0;
    double ds = 0.0;
    for (int i = 0; i < 3; i++) {
      dt = fmax(dt, hypot(strobe_t->data[i].re - direct->data[i].re,
                          strobe_t->data[i].im - direct->data[i].im));
      ds = fmax(ds, hypot(strobe_s->data[i].re - direct->data[i].re,
                          strobe_s->data[i].im - direct->data[i].im));
    }

    check_close(dt, 0.0, 1e-9,
                "time-domain stroboscopic vs direct (7 periods)");
    check_close(ds, 0.0, 1e-7, "Sambe stroboscopic vs direct (7 periods)");
    check_close(cvector_norm(strobe_t), 1.0, 1e-10, "stroboscopic norm");
  } else {
    check_true(0, "stroboscopic states allocated");
  }

  cvector_t *zero = floquet_stroboscopic_state(rt, psi0, 0);
  if (zero) {
    check_close(hypot(zero->data[1].re - psi0->data[1].re,
                      zero->data[1].im - psi0->data[1].im),
                0.0, 1e-10, "n_periods = 0 returns psi0");

    cvector_free(zero);
  }

  cvector_free(strobe_t);
  cvector_free(strobe_s);
  cvector_free(direct);
  cvector_free(psi0);
  floquet_result_free(rs);

  // Truncation weight falls as K grows
  floquet_result_t *k2 = floquet_solve_sambe(3, h, 1, p.om, 2);
  floquet_result_t *k5 = floquet_solve_sambe(3, h, 1, p.om, 5);
  if (k2 && k5) {
    printf("  truncation weight: K=2 -> %.3e, K=5 -> %.3e\n",
           k2->truncation_weight, k5->truncation_weight);

    check_true(k5->truncation_weight < k2->truncation_weight,
               "edge weight decreases with K");
  } else {
    check_true(0, "small-K Sambe solves returned results");
  }

  floquet_result_free(k2);
  floquet_result_free(k5);
  floquet_result_free(rt);
  free_harmonics(h, 3);
}

static void test_static_and_degenerate(void) {
  printf("  === Static and degenerate limits ===\n");

  static_params_t sp = {0.9, -0.2};
  floquet_result_t *r = floquet_solve_time(2, static_H, &sp, 1.0, 2000);
  if (r) {
    // 0.9 folds to -0.1 (mod 1); -0.2 is already inside [-0.5, 0.5)
    check_close(r->quasienergies[0], -0.2, 1e-10, "static: -0.2");
    check_close(r->quasienergies[1], -0.1, 1e-10, "static: 0.9 folded to -0.1");

    floquet_result_free(r);
  } else {
    check_true(0, "static time solve returned a result");
  }

  cmatrix_t *h0 = cmatrix_alloc(2, 2);
  put(h0, 0, 0, sp.eps0, 0.0);
  put(h0, 1, 1, sp.eps1, 0.0);

  cmatrix_t *hs[1] = {h0};
  floquet_result_t *rs = floquet_solve_sambe(2, hs, 0, 1.0, 3);
  if (rs) {
    check_close(rs->quasienergies[0], -0.2, 1e-10, "Sambe static: -0.2");
    check_close(rs->quasienergies[1], -0.1, 1e-10, "Sambe static: -0.1");

    floquet_result_free(rs);
  } else {
    check_true(0, "static Sambe solve returned a result");
  }

  cmatrix_free(h0);

  // H = 0: fully degenerate, U = I
  floquet_result_t *rz = floquet_solve_time(3, zero_H, NULL, 2.0, 16);
  if (rz) {
    for (int a = 0; a < 3; a++) {
      char label[48];
      snprintf(label, sizeof label, "H=0 quasi-energy %d", a);

      check_close(rz->quasienergies[a], 0.0, 1e-12, label);
    }

    check_close(orthonormality_error(rz->modes), 0.0, 1e-10,
                "H=0 modes orthonormal despite full degeneracy");

    floquet_result_free(rz);
  } else {
    check_true(0, "H=0 solve returned a result");
  }
}

static void test_invalid(void) {
  printf("  === Invalid input ===\n");

  circ_params_t p = {1.0, 0.3, 0.85};
  check_true(floquet_solve_time(0, circular_H, &p, 1.0, 100) == NULL, "n = 0");
  check_true(floquet_solve_time(2, NULL, &p, 1.0, 100) == NULL,
             "NULL callback");
  check_true(floquet_solve_time(2, circular_H, &p, 0.0, 100) == NULL,
             "\\omega = 0");
  check_true(floquet_solve_time(2, circular_H, &p, -1.0, 100) == NULL,
             "\\omega < 0");
  check_true(floquet_solve_time(2, circular_H, &p, NAN, 100) == NULL,
             "\\omega = NaN");
  check_true(floquet_solve_time(2, circular_H, &p, 1.0, 3) == NULL,
             "steps < 4");
  check_true(floquet_solve_time(100000, circular_H, &p, 1.0, 100) == NULL,
             "n above supported limit");

  cmatrix_t *h[3];
  circular_harmonics(&p, h);
  check_true(floquet_solve_sambe(2, NULL, 1, 1.0, 4) == NULL, "NULL harmonics");
  check_true(floquet_solve_sambe(2, h, -1, 1.0, 4) == NULL, "M < 0");
  check_true(floquet_solve_sambe(2, h, 1, 1.0, 1) == NULL, "K < M + 1");
  check_true(floquet_solve_sambe(2, h, 1, 0.0, 4) == NULL, "Sambe \\omega = 0");
  check_true(floquet_solve_sambe(2, h, 1, 1.0, 100000) == NULL,
             "oversized extended space rejected");
  check_true(floquet_solve_sambe(3, h, 1, 1.0, 4) == NULL,
             "dimension mismatch with harmonics");

  // Non-Hermitian drive: H_{-1} != H_{+1}^\dagger
  put(h[0], 0, 1, 0.5, 0.0);
  put(h[2], 1, 0, 0.9, 0.0);
  check_true(floquet_solve_sambe(2, h, 1, 1.0, 4) == NULL,
             "non-Hermitian harmonics rejected");

  free_harmonics(h, 3);

  floquet_result_t *r = floquet_solve_time(2, circular_H, &p, p.om, 200);
  cvector_t *wrong = cvector_alloc(3);

  check_true(floquet_stroboscopic_state(NULL, wrong, 1) == NULL,
             "strobe: NULL result");
  check_true(floquet_stroboscopic_state(r, NULL, 1) == NULL,
             "strobe: NULL psi");
  check_true(floquet_stroboscopic_state(r, wrong, 1) == NULL,
             "strobe: wrong length");

  cvector_t *ok = cvector_alloc(2);
  check_true(floquet_stroboscopic_state(r, ok, -1) == NULL,
             "strobe: negative period count");

  cvector_free(ok);
  cvector_free(wrong);
  floquet_result_free(r);
  floquet_result_free(NULL);

  check_true(isnan(floquet_circular_rabi_quasienergy(1.0, 0.3, 0.0, 1)),
             "Rabi helper: \\omega = 0 -> NaN");
  check_true(isnan(floquet_circular_rabi_quasienergy(1.0, 0.3, 1.0, 0)),
             "Rabi helper: bad sign -> NaN");
}

int main(void) {
  test_bessel();
  test_circular_rabi();
  test_rk4_order();
  test_cdt();
  test_bloch_siegert();
  test_spin1();
  test_static_and_degenerate();
  test_invalid();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll Floquet checks passed\n");
  return 0;
}
