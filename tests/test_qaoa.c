/*
 * Test: Quantum Approximate Optimization Algorithm (QAOA).
 *
 * 1. Cost tables: MaxCut on K3/K4/ring against known cut values; weighted
 *    edges; Ising table against a hand-evaluated spin configuration
 * 2. Trivial limits: gamma = 0 or beta = 0 leaves uniform superposition,
 *    so <C> equals mean of cost table; states stay normalized.
 * 3. Exact closed forms: ring of n vertices at p = 1 has
 *    <C> = n (1/2 + \sin(4\beta) sin(2\gamma) / 4); general-graph p = 1
 *    formula (triangles included) matches state-vector simulator on
 *    random graphs
 * 4. Cross-check against dense Pauli-Hamiltonian path in vqe.c: build
 *    C = \sum w (I - Z_u Z_v) / 2 as a matrix and compare <\psi|C|\psi>
 * 5. Periodicity: integer costs make <C> periodic in \gamma (2\pi) and \beta
 *    (\pi)
 * 6. Optimizer: on a ring optimal depth-p expectation is exactly
 *    (2p + 1)/(2p + 2) of edges (Farhi et al.), valid for n > 2p + 1;
 *    optimum probability grows with p; Ising minimization works. On
 *    Petersen graph p = 1 optimum is |E| (1/2 + 1/(3\sqrt(3)))
 * 7. Invalid-input handling
 */

#include "../core/matrix.h"
#include "../core/random.h"
#include "../physics/qaoa.h"
#include "../physics/vqe.h"
#include "vector.h"
#include <math.h>
#include <stdint.h>
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
  printf("  %s: got=%.10f expected=%.10f err=%.2e\n", label, got, expected,
         err);
  if (err > tol) {
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

static qaoa_problem_t *make_ring(int n) {
  int *u = malloc((size_t)n * sizeof *u);
  int *v = malloc((size_t)n * sizeof *v);
  qaoa_problem_t *pr = NULL;
  if (u && v) {
    for (int i = 0; i < n; i++) {
      u[i] = i;
      v[i] = (i + 1) % n;
    }

    pr = qaoa_problem_maxcut(n, u, v, NULL, n);
  }

  free(u);
  free(v);

  return pr;
}

static void test_cost_tables(void) {
  printf("  === Cost tables and exact optima ===\n");

  const int tu[] = {0, 1, 0};
  const int tv[] = {1, 2, 2};
  qaoa_problem_t *k3 = qaoa_problem_maxcut(3, tu, tv, NULL, 3);
  check_true(k3 != NULL, "K3 built");
  int best = -1;
  check_close(qaoa_problem_best_cost(k3, &best), 2.0, 1e-12, "K3 max cut");

  check_true(best >= 1 && best <= 6, "K3 best bitstring is a real cut");

  check_close(qaoa_problem_worst_cost(k3), 0.0, 1e-12, "K3 min cut");
  check_close(qaoa_approximation_ratio(k3, 1.0), 0.5, 1e-12,
              "K3 ratio at <C>=1");
  check_true(qaoa_problem_qubits(k3) == 3, "K3 qubit count");

  qaoa_problem_free(k3);

  int ku[6], kv[6], m = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = i + 1; j < 4; j++) {
      ku[m] = i;
      kv[m] = j;
      m++;
    }
  }

  qaoa_problem_t *k4 = qaoa_problem_maxcut(4, ku, kv, NULL, 6);
  check_close(qaoa_problem_best_cost(k4, NULL), 4.0, 1e-12, "K4 max cut");

  qaoa_problem_free(k4);

  qaoa_problem_t *ring = make_ring(8);
  check_close(qaoa_problem_best_cost(ring, NULL), 8.0, 1e-12,
              "ring(8) max cut = all edges (bipartite)");

  qaoa_problem_free(ring);

  // Weighted: path 0-1-2 with weights 2 and 3; cutting both = 5
  const int wu[] = {0, 1};
  const int wv[] = {1, 2};
  const double ww[] = {2.0, 3.0};
  qaoa_problem_t *wp = qaoa_problem_maxcut(3, wu, wv, ww, 2);
  check_close(qaoa_problem_best_cost(wp, NULL), 5.0, 1e-12,
              "weighted path max cut");

  const double *tab = qaoa_problem_cost_table(wp);

  // index 0b010 = qubits (0,1,2) = (0,1,0): both edges cut
  check_close(tab[2], 5.0, 1e-12, "table[010] (qubit 0 = MSB)");

  // index 0b100 = (1,0,0): only edge (0,1) cut, weight 2
  check_close(tab[4], 2.0, 1e-12, "table[100]");

  qaoa_problem_free(wp);

  // Ising: C = 0.5 + 1.0*z0 - 2.0*z1 + 0.7*z0 z1
  const double h[] = {1.0, -2.0};
  const int ci[] = {0};
  const int cj[] = {1};
  const double J[] = {0.7};
  qaoa_problem_t *is = qaoa_problem_ising(2, h, ci, cj, J, 1, 0.5);

  check_true(is != NULL, "Ising built");

  tab = qaoa_problem_cost_table(is);

  // index 0b01 => z0=+1, z1=-1: 0.5 + 1 + 2 - 0.7 = 2.8
  check_close(tab[1], 2.8, 1e-12, "Ising table[01]");

  // index 0b10 => z0=-1, z1=+1: 0.5 - 1 - 2 - 0.7 = -3.2 (minimum)
  int ib = -1;
  check_close(qaoa_problem_best_cost(is, &ib), -3.2, 1e-12,
              "Ising minimum energy");
  check_true(ib == 2, "Ising argmin is 0b10");

  qaoa_problem_free(is);
}

static void test_trivial_limits(void) {
  printf("  === Trivial limits and normalization ===\n");

  qaoa_problem_t *ring = make_ring(6);
  const double *tab = qaoa_problem_cost_table(ring);

  double mean = 0.0;
  for (int k = 0; k < 64; k++) {
    mean += tab[k];
  }

  mean /= 64.0;
  check_close(mean, 3.0, 1e-12, "mean cut of ring(6) = |E|/2");

  double g = 0.9;
  double b = 0.0;

  check_close(qaoa_expectation(ring, 1, &g, &b), mean, 1e-12, "beta = 0");
  g = 0.0;
  b = 0.4;

  check_close(qaoa_expectation(ring, 1, &g, &b), mean, 1e-12, "gamma = 0");

  const double gs[3] = {0.3, -1.1, 0.8};
  const double bs[3] = {0.5, 0.2, -0.7};
  cvector_t *psi = qaoa_prepare_state(ring, 3, gs, bs);
  check_true(psi != NULL, "p=3 state prepared");
  if (psi) {
    check_close(cvector_norm(psi), 1.0, 1e-12, "state norm (p=3)");

    cvector_free(psi);
  }

  qaoa_problem_free(ring);
}

static void test_closed_forms(void) {
  printf("  === Closed-form p = 1 expectations ===\n");

  qaoa_problem_t *ring = make_ring(8);
  double g = 0.7, b = 0.3;
  double expected = 8.0 * (0.5 + 0.25 * sin(4.0 * b) * sin(2.0 * g));
  check_close(qaoa_expectation(ring, 1, &g, &b), expected, 1e-11,
              "ring(8) simulator vs closed form");
  check_close(qaoa_maxcut_p1_expectation(ring, g, b), expected, 1e-11,
              "ring(8) Wang formula vs closed form");

  qaoa_problem_free(ring);

  rng_state_t rng;
  rng_seed(&rng, 12345);
  for (int trial = 0; trial < 4; trial++) {
    int n = 7;
    int u[64];
    int v[64];
    int m = 0;

    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        if (rng_uniform(&rng) < 0.5) {
          u[m] = i;
          v[m] = j;
          m++;
        }
      }
    }

    qaoa_problem_t *pr = qaoa_problem_maxcut(n, u, v, NULL, m);
    double gg = rng_uniform_range(&rng, -3.0, 3.0);
    double bb = rng_uniform_range(&rng, -1.5, 1.5);
    char label[64];
    snprintf(label, sizeof label, "random graph %d (%d edges)", trial, m);
    check_close(qaoa_maxcut_p1_expectation(pr, gg, bb),
                qaoa_expectation(pr, 1, &gg, &bb), 1e-11, label);

    qaoa_problem_free(pr);
  }

  // Formula does not apply to weighted problems
  const int wu[] = {0, 1};
  const int wv[] = {1, 2};
  const double ww[] = {1.0, 2.0};
  qaoa_problem_t *wp = qaoa_problem_maxcut(3, wu, wv, ww, 2);
  check_true(isnan(qaoa_maxcut_p1_expectation(wp, 0.3, 0.2)),
             "weighted problem -> NaN");
  qaoa_problem_free(wp);

  const int is_ci[] = {0};
  const int is_cj[] = {1};
  const double is_J[] = {1.0};
  qaoa_problem_t *is = qaoa_problem_ising(2, NULL, is_ci, is_cj, is_J, 1, 0.0);
  check_true(isnan(qaoa_maxcut_p1_expectation(is, 0.3, 0.2)),
             "Ising problem -> NaN");

  qaoa_problem_free(is);
}

static void test_dense_crosscheck(void) {
  printf("  === Cross-check against dense Pauli Hamiltonian (vqe.c) ===\n");
  // Triangle with a pendant vertex: 4 qubits, weights all 1
  int n = 4;
  const int u[] = {0, 1, 0, 2};
  const int v[] = {1, 2, 2, 3};
  int m = 4;
  qaoa_problem_t *pr = qaoa_problem_maxcut(n, u, v, NULL, m);

  // C = \sum (I - Z_u Z_v)/2  =>  (m/2) I - (1/2) \sum Z_u Z_v
  char strs[5][8];
  const char *ptrs[5];
  double coeff[5];
  snprintf(strs[0], sizeof strs[0], "IIII");
  coeff[0] = 0.5 * m;
  for (int e = 0; e < m; e++) {
    memset(strs[e + 1], 'I', (size_t)n);
    strs[e + 1][u[e]] = 'Z';
    strs[e + 1][v[e]] = 'Z';
    strs[e + 1][n] = '\0';
    coeff[e + 1] = -0.5;
  }

  for (int i = 0; i < 5; i++) {
    ptrs[i] = strs[i];
  }

  cmatrix_t *H = vqe_build_pauli_hamiltonian(n, ptrs, coeff, 5);

  check_true(H != NULL, "dense C built");

  const double gs[2] = {0.42, -0.9};
  const double bs[2] = {0.31, 0.66};
  cvector_t *psi = qaoa_prepare_state(pr, 2, gs, bs);
  if (H && psi) {
    check_close(vqe_expectation(psi, H), qaoa_expectation(pr, 2, gs, bs), 1e-11,
                "dense <C> vs diagonal <C> (p=2)");
  }

  cvector_free(psi);
  cmatrix_free(H);
  qaoa_problem_free(pr);
}

static void test_periodicity(void) {
  printf("  === Periodicity (integer costs) ===\n");

  qaoa_problem_t *ring = make_ring(8);
  double g0 = 0.7, b0 = 0.3;
  double g1 = g0 + 2.0 * M_PI;
  double b1 = b0 + M_PI;

  check_close(qaoa_expectation(ring, 1, &g1, &b1),
              qaoa_expectation(ring, 1, &g0, &b0), 1e-9,
              "<C>(\\gamma + 2 \\pi, \\beta + \\pi) = <C>(\\gamma, \\beta)");

  qaoa_problem_free(ring);
}

static void test_optimizer(void) {
  printf("  === Optimizer ===\n");

  // Farhi et al.: ring, p layers, optimal <C>/edges = (2p+1)/(2p+2) for
  // n > 2p + 1. n = 8 covers p = 1, 2, 3
  qaoa_problem_t *ring = make_ring(8);
  double prev_prob = 0.0;
  for (int p = 1; p <= 3; p++) {
    qaoa_result_t r = qaoa_run(ring, p, 8, 6, 0.4, 2026);

    check_true(r.gamma != NULL && r.beta != NULL, "result allocated");
    char label[80];

    snprintf(label, sizeof label, "ring(8) p=%d optimal <C>/edges", p);
    check_close(r.expectation / 8.0, (2.0 * p + 1.0) / (2.0 * p + 2.0), 1e-8,
                label);

    snprintf(label, sizeof label, "ring(8) p=%d ratio == <C>/8", p);
    check_close(r.approx_ratio, r.expectation / 8.0, 1e-12, label);

    snprintf(label, sizeof label, "ring(8) p=%d optimum prob increased", p);
    check_true(r.optimum_prob > prev_prob, label);

    prev_prob = r.optimum_prob;
    if (p == 3) {
      // best bitstring must be one of two perfect alternating cuts
      check_true(r.best_bitstring == 0x55 || r.best_bitstring == 0xAA,
                 "most probable bitstring is an alternating cut");
    }

    qaoa_result_free(&r);
  }

  qaoa_problem_free(ring);

  // Deterministic: same seed -> identical result
  ring = make_ring(6);
  qaoa_result_t a = qaoa_run(ring, 2, 4, 4, 0.4, 7);
  qaoa_result_t b = qaoa_run(ring, 2, 4, 4, 0.4, 7);
  check_close(a.expectation, b.expectation, 0.0, "same seed reproducible");
  qaoa_result_free(&a);
  qaoa_result_free(&b);
  qaoa_problem_free(ring);

  // Ising minimization: 1D antiferromagnetic chain with a field; expectation
  // must drop below uniform-superposition mean and approach ground energy as
  // depth grows
  int n = 6;
  int ci[5];
  int cj[5];
  double J[5];
  for (int i = 0; i < 5; i++) {
    ci[i] = i;
    cj[i] = i + 1;
    J[i] = 1.0;
  }

  const double h[6] = {0.3, -0.2, 0.1, 0.0, -0.1, 0.2};
  qaoa_problem_t *is = qaoa_problem_ising(n, h, ci, cj, J, 5, 0.0);
  const double *tab = qaoa_problem_cost_table(is);
  double mean = 0.0;
  for (int k = 0; k < 64; k++) {
    mean += tab[k];
  }
  mean /= 64.0;

  double emin = qaoa_problem_best_cost(is, NULL);
  qaoa_result_t r1 = qaoa_run(is, 1, 8, 6, 0.4, 99);
  qaoa_result_t r3 = qaoa_run(is, 3, 8, 6, 0.4, 99);

  printf("  Ising: mean=%.6f  E0=%.6f  <C>_p1=%.6f  <C>_p3=%.6f\n", mean, emin,
         r1.expectation, r3.expectation);
  check_true(r1.expectation < mean - 1e-3, "p=1 beats uniform superposition");
  check_true(r3.expectation < r1.expectation + 1e-9,
             "p=3 no worse than p=1 (minimization)");
  check_true(r3.expectation >= emin - 1e-9, "variational bound <C> >= E0");
  check_true(r3.approx_ratio > r1.approx_ratio - 1e-9 &&
                 r3.approx_ratio <= 1.0 + 1e-12,
             "ratio in range and non-decreasing with depth");

  qaoa_result_free(&r1);
  qaoa_result_free(&r3);
  qaoa_problem_free(is);
}

static void test_petersen(void) {
  printf("  === Petersen graph (triangle-free, 3-regular) ===\n");
  const int u[15] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 7, 9, 6, 8};
  const int v[15] = {1, 2, 3, 4, 0, 5, 6, 7, 8, 9, 7, 9, 6, 8, 5};
  qaoa_problem_t *pr = qaoa_problem_maxcut(10, u, v, NULL, 15);
  check_true(pr != NULL, "Petersen built");
  if (!pr) {
    return;
  }
  check_close(qaoa_problem_best_cost(pr, NULL), 12.0, 1e-12,
              "Petersen exact max cut");

  // d = e = 2, f = 0 on every edge, so <C_uv> = 1/2 + \sin(4 \beta)
  // \sin(\gamma)
  // \cos^2(\gamma) / 2, maximized at beta = pi/8 and gamma maximizing
  // \sin(\gamma) \cos^2(\gamma) (= 2/(3 sqrt 3)): |E| (1/2 + 1/(3\sqrt(3)))
  double p1_exact = 15.0 * (0.5 + 1.0 / (3.0 * sqrt(3.0)));
  qaoa_result_t r = qaoa_run(pr, 1, 4, 6, 0.4, 2026);
  check_close(r.expectation, p1_exact, 1e-8,
              "Petersen p=1 optimum vs closed form");
  check_close(qaoa_maxcut_p1_expectation(pr, r.gamma[0], r.beta[0]),
              r.expectation, 1e-10, "Wang formula at optimized angles");
  qaoa_result_free(&r);
  qaoa_result_t r2 = qaoa_run(pr, 2, 4, 6, 0.4, 2026);
  check_true(r2.expectation > p1_exact + 1e-3, "p=2 improves on p=1");
  check_true(r2.expectation <= 12.0 + 1e-9, "<C> never exceeds max cut");
  qaoa_result_free(&r2);
  qaoa_problem_free(pr);
}

static void test_invalid_input(void) {
  printf("  === Invalid input ===\n");

  const int u[] = {0, 1};
  const int v[] = {1, 2};
  const int selfu[] = {0};
  const int selfv[] = {0};
  const int oobu[] = {0};
  const int oobv[] = {5};
  const int okv[] = {1};
  const double nanw[] = {NAN};

  check_true(qaoa_problem_maxcut(0, u, v, NULL, 2) == NULL, "n = 0");
  check_true(qaoa_problem_maxcut(QAOA_MAX_QUBITS + 1, u, v, NULL, 2) == NULL,
             "n too large");
  check_true(qaoa_problem_maxcut(3, NULL, v, NULL, 2) == NULL, "NULL u");
  check_true(qaoa_problem_maxcut(3, u, v, NULL, 0) == NULL, "zero edges");
  check_true(qaoa_problem_maxcut(3, selfu, selfv, NULL, 1) == NULL,
             "self-loop");
  check_true(qaoa_problem_maxcut(3, oobu, oobv, NULL, 1) == NULL,
             "out-of-range endpoint");
  check_true(qaoa_problem_maxcut(3, selfu, okv, nanw, 1) == NULL,
             "non-finite weight");
  check_true(qaoa_problem_ising(2, NULL, NULL, NULL, NULL, 1, 0.0) == NULL,
             "Ising: couplings without arrays");
  check_true(qaoa_problem_ising(2, NULL, NULL, NULL, NULL, -1, 0.0) == NULL,
             "Ising: negative coupling count");
  check_true(qaoa_problem_ising(2, NULL, NULL, NULL, NULL, 0, NAN) == NULL,
             "Ising: NaN offset");
  qaoa_problem_t *fields_only =
      qaoa_problem_ising(2, (double[]){1.0, -1.0}, NULL, NULL, NULL, 0, 0.0);
  check_true(fields_only != NULL, "Ising: fields only is valid");
  qaoa_problem_free(fields_only);

  qaoa_problem_t *ring = make_ring(4);
  double g = 0.1, b = 0.2;
  check_true(qaoa_prepare_state(NULL, 1, &g, &b) == NULL, "NULL problem");
  check_true(qaoa_prepare_state(ring, 0, &g, &b) == NULL, "p = 0");
  check_true(qaoa_prepare_state(ring, 1, NULL, &b) == NULL, "NULL gamma");
  double bad = NAN;
  check_true(qaoa_prepare_state(ring, 1, &bad, &b) == NULL, "NaN gamma");
  check_close(qaoa_expectation(NULL, 1, &g, &b), 0.0, 0.0,
              "expectation of NULL problem is 0");
  qaoa_result_t r = qaoa_run(ring, 0, 1, 1, 0.3, 1);
  check_true(r.gamma == NULL && r.beta == NULL, "run p = 0 rejected");
  r = qaoa_run(ring, 1, 0, 1, 0.3, 1);
  check_true(r.gamma == NULL, "run n_restarts = 0 rejected");
  r = qaoa_run(ring, 1, 1, 1, -0.3, 1);
  check_true(r.gamma == NULL, "run negative window rejected");
  qaoa_result_free(&r);
  qaoa_problem_free(ring);
  qaoa_problem_free(NULL);
  check_true(isnan(qaoa_maxcut_p1_expectation(NULL, 0.1, 0.1)),
             "p1 formula on NULL -> NaN");
  check_close(qaoa_approximation_ratio(NULL, 1.0), 0.0, 0.0,
              "ratio of NULL problem is 0");
}

int main(void) {
  printf("  >  Test QAOA :\n");

  test_cost_tables();
  test_trivial_limits();
  test_closed_forms();
  test_dense_crosscheck();
  test_periodicity();
  test_optimizer();
  test_petersen();
  test_invalid_input();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll QAOA checks passed\n");
  return 0;
}
