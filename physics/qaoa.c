#include "qaoa.h"
#include "../core/complex.h"
#include "../core/random.h"
#include "qubits.h"
#include "variational.h"
#include "vector.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct qaoa_problem {
  int n;        // number of qubits
  int dim;      // 2^n
  int maximize; // 1: maximize C, 0: minimize C
  double *cost; // length dim
  // MaxCut bookkeeping, used only by qaoa_maxcut_p1_expectation
  int is_maxcut;
  int unit_simple; // MaxCut with all weights == 1 and no parallel edges
  int n_edges;
  int *eu;
  int *ev;
};

static int qaoa_bit(int n, int k, int qubit) {
  return (k >> (n - 1 - qubit)) & 1;
}

static int qaoa_popcount(uint32_t x) {
  int c = 0;
  while (x) {
    x &= x - 1U;
    c++;
  }

  return c;
}

static qaoa_problem_t *qaoa_problem_alloc(int n) {
  qaoa_problem_t *pr = calloc(1, sizeof *pr);
  if (!pr) {
    return NULL;
  }

  pr->n = n;
  pr->dim = 1 << n;
  pr->cost = calloc((size_t)pr->dim, sizeof *pr->cost);
  if (!pr->cost) {
    free(pr);

    return NULL;
  }

  return pr;
}

qaoa_problem_t *qaoa_problem_maxcut(int n_qubits, const int *u, const int *v,
                                    const double *w, int n_edges) {
  if (n_qubits < 1 || n_qubits > QAOA_MAX_QUBITS || !u || !v || n_edges < 1) {
    return NULL;
  }

  for (int e = 0; e < n_edges; e++) {
    if (u[e] < 0 || u[e] >= n_qubits || v[e] < 0 || v[e] >= n_qubits ||
        u[e] == v[e] || (w && !isfinite(w[e]))) {
      return NULL;
    }
  }

  qaoa_problem_t *pr = qaoa_problem_alloc(n_qubits);
  if (!pr) {
    return NULL;
  }

  pr->maximize = 1;
  pr->is_maxcut = 1;
  pr->n_edges = n_edges;
  pr->eu = malloc((size_t)n_edges * sizeof *pr->eu);
  pr->ev = malloc((size_t)n_edges * sizeof *pr->ev);
  if (!pr->eu || !pr->ev) {
    qaoa_problem_free(pr);

    return NULL;
  }

  uint32_t adj[QAOA_MAX_QUBITS] = {0};
  int unit_simple = 1;
  for (int e = 0; e < n_edges; e++) {
    double we = w ? w[e] : 1.0;

    pr->eu[e] = u[e];
    pr->ev[e] = v[e];
    if (we != 1.0 || (adj[u[e]] & (1U << v[e]))) {
      unit_simple = 0;
    }

    adj[u[e]] |= 1U << v[e];
    adj[v[e]] |= 1U << u[e];

    for (int k = 0; k < pr->dim; k++) {
      if (qaoa_bit(n_qubits, k, u[e]) != qaoa_bit(n_qubits, k, v[e])) {
        pr->cost[k] += we;
      }
    }
  }

  pr->unit_simple = unit_simple;

  return pr;
}

qaoa_problem_t *qaoa_problem_ising(int n_qubits, const double *h, const int *ci,
                                   const int *cj, const double *J,
                                   int n_couplings, double offset) {
  if (n_qubits < 1 || n_qubits > QAOA_MAX_QUBITS || n_couplings < 0 ||
      !isfinite(offset) || (n_couplings > 0 && (!ci || !cj || !J))) {
    return NULL;
  }

  if (h) {
    for (int i = 0; i < n_qubits; i++) {
      if (!isfinite(h[i])) {
        return NULL;
      }
    }
  }

  for (int c = 0; c < n_couplings; c++) {
    if (ci[c] < 0 || ci[c] >= n_qubits || cj[c] < 0 || cj[c] >= n_qubits ||
        ci[c] == cj[c] || !isfinite(J[c])) {
      return NULL;
    }
  }

  qaoa_problem_t *pr = qaoa_problem_alloc(n_qubits);
  if (!pr) {
    return NULL;
  }

  pr->maximize = 0;

  for (int k = 0; k < pr->dim; k++) {
    pr->cost[k] = offset;
  }

  if (h) {
    for (int i = 0; i < n_qubits; i++) {
      for (int k = 0; k < pr->dim; k++) {
        pr->cost[k] += h[i] * (1.0 - 2.0 * qaoa_bit(n_qubits, k, i));
      }
    }
  }

  for (int c = 0; c < n_couplings; c++) {
    for (int k = 0; k < pr->dim; k++) {
      double zi = 1.0 - 2.0 * qaoa_bit(n_qubits, k, ci[c]);
      double zj = 1.0 - 2.0 * qaoa_bit(n_qubits, k, cj[c]);

      pr->cost[k] += J[c] * zi * zj;
    }
  }

  return pr;
}

void qaoa_problem_free(qaoa_problem_t *problem) {
  if (!problem) {
    return;
  }

  free(problem->cost);
  free(problem->eu);
  free(problem->ev);
  free(problem);
}

int qaoa_problem_qubits(const qaoa_problem_t *problem) {
  return problem ? problem->n : 0;
}

const double *qaoa_problem_cost_table(const qaoa_problem_t *problem) {
  return problem ? problem->cost : NULL;
}

double qaoa_problem_best_cost(const qaoa_problem_t *problem, int *best_index) {
  if (!problem) {
    return 0.0;
  }

  int best = 0;
  for (int k = 1; k < problem->dim; k++) {
    if (problem->maximize ? problem->cost[k] > problem->cost[best]
                          : problem->cost[k] < problem->cost[best]) {
      best = k;
    }
  }

  if (best_index) {
    *best_index = best;
  }

  return problem->cost[best];
}

double qaoa_problem_worst_cost(const qaoa_problem_t *problem) {
  if (!problem) {
    return 0.0;
  }

  double worst = problem->cost[0];
  for (int k = 1; k < problem->dim; k++) {
    if (problem->maximize ? problem->cost[k] < worst
                          : problem->cost[k] > worst) {
      worst = problem->cost[k];
    }
  }

  return worst;
}

static int qaoa_angles_ok(int p, const double *gamma, const double *beta) {
  if (p < 1 || !gamma || !beta) {
    return 0;
  }

  for (int l = 0; l < p; l++) {
    if (!isfinite(gamma[l]) || !isfinite(beta[l])) {
      return 0;
    }
  }

  return 1;
}

cvector_t *qaoa_prepare_state(const qaoa_problem_t *problem, int p,
                              const double *gamma, const double *beta) {
  if (!problem || !qaoa_angles_ok(p, gamma, beta)) {
    return NULL;
  }

  cvector_t *psi = cvector_alloc(problem->dim);
  if (!psi) {
    return NULL;
  }

  double amp = 1.0 / sqrt((double)problem->dim);
  for (int k = 0; k < problem->dim; k++) {
    psi->data[k] = c_real(amp);
  }

  complex_t rx[4];
  for (int l = 0; l < p; l++) {
    // Phase separator e^{-i \gamma C}: diagonal in computational basis
    for (int k = 0; k < problem->dim; k++) {
      double phi = -gamma[l] * problem->cost[k];

      psi->data[k] = c_mul(psi->data[k], c_new(cos(phi), sin(phi)));
    }

    // Mixer e^{-i \beta \sum X_i} = \prod_i RX_i(2 \beta)
    rx_gate(2.0 * beta[l], rx);
    for (int q = 0; q < problem->n; q++) {
      qstate_apply_gate1(psi, problem->n, q, rx);
    }
  }

  return psi;
}

static double qaoa_expect_state(const qaoa_problem_t *problem,
                                const cvector_t *psi) {
  double e = 0.0;
  for (int k = 0; k < problem->dim; k++) {
    e += c_abs2(psi->data[k]) * problem->cost[k];
  }

  return e;
}

double qaoa_expectation(const qaoa_problem_t *problem, int p,
                        const double *gamma, const double *beta) {
  cvector_t *psi = qaoa_prepare_state(problem, p, gamma, beta);
  if (!psi) {
    return 0.0;
  }

  double e = qaoa_expect_state(problem, psi);

  cvector_free(psi);

  return e;
}

double qaoa_maxcut_p1_expectation(const qaoa_problem_t *problem, double gamma,
                                  double beta) {
  if (!problem || !problem->is_maxcut || !problem->unit_simple ||
      !isfinite(gamma) || !isfinite(beta)) {
    return NAN;
  }

  uint32_t adj[QAOA_MAX_QUBITS] = {0};
  for (int e = 0; e < problem->n_edges; e++) {
    adj[problem->eu[e]] |= 1U << problem->ev[e];
    adj[problem->ev[e]] |= 1U << problem->eu[e];
  }

  double cg = cos(gamma);
  double s4b = sin(4.0 * beta);
  double s2b = sin(2.0 * beta);
  double sg = sin(gamma);
  double c2g = cos(2.0 * gamma);

  double total = 0.0;
  for (int e = 0; e < problem->n_edges; e++) {
    int u = problem->eu[e];
    int v = problem->ev[e];
    int d = qaoa_popcount(adj[u]) - 1;
    int ee = qaoa_popcount(adj[v]) - 1;
    int f = qaoa_popcount(adj[u] & adj[v]);

    total += 0.5 + 0.25 * s4b * sg * (pow(cg, d) + pow(cg, ee)) -
             0.25 * s2b * s2b * pow(cg, d + ee - 2 * f) * (1.0 - pow(c2g, f));
  }

  return total;
}

double qaoa_approximation_ratio(const qaoa_problem_t *problem,
                                double expectation) {
  if (!problem) {
    return 0.0;
  }

  double best = qaoa_problem_best_cost(problem, NULL);
  double worst = qaoa_problem_worst_cost(problem);
  if (best == worst) {
    return 1.0;
  }

  return (expectation - worst) / (best - worst);
}

typedef struct {
  const qaoa_problem_t *problem;
  int p;
  double *x; // 2p angles: gamma[0..p-1], then beta[0..p-1]
  int idx;   // coordinate being varied
} qaoa_coord_closure_t;

// Objective in "minimize" form: -<C> for maximization problems, <C> otherwise
static double qaoa_objective(const qaoa_problem_t *problem, int p,
                             const double *x) {
  double e = qaoa_expectation(problem, p, x, x + p);

  return problem->maximize ? -e : e;
}

static double qaoa_coord_eval(double val, void *params) {
  qaoa_coord_closure_t *cl = params;
  double saved = cl->x[cl->idx];
  cl->x[cl->idx] = val;

  double f = qaoa_objective(cl->problem, cl->p, cl->x);
  cl->x[cl->idx] = saved;

  return f;
}

// Standard deviation of cost table: natural energy scale for setting initial
// phase-separator angle
static double qaoa_cost_sigma(const qaoa_problem_t *problem) {
  double mean = 0.0;
  for (int k = 0; k < problem->dim; k++) {
    mean += problem->cost[k];
  }

  mean /= (double)problem->dim;
  double var = 0.0;
  for (int k = 0; k < problem->dim; k++) {
    double d = problem->cost[k] - mean;
    var += d * d;
  }

  var /= (double)problem->dim;
  double sigma = sqrt(var);

  return sigma > 1e-12 ? sigma : 1.0;
}

// Central-difference gradient of objective (2p angles, 4p evaluations)
static void qaoa_gradient(const qaoa_problem_t *problem, int p, double *x,
                          double *g) {
  const double h = 1e-6;
  for (int i = 0; i < 2 * p; i++) {
    double xi = x[i];
    x[i] = xi + h;

    double fp = qaoa_objective(problem, p, x);
    x[i] = xi - h;

    double fm = qaoa_objective(problem, p, x);
    x[i] = xi;
    g[i] = (fp - fm) / (2.0 * h);
  }
}

// BFGS refinement (backtracking Armijo line search) of 2p angles in x.
// Coordinate descent stalls in narrow, strongly correlated gamma-beta
// valleys of QAOA landscape at p >= 2; a quasi-Newton polish reaches
// local optimum to ~1e-10. Returns final objective. On any allocation
// failure x is left untouched at its last accepted value
static double qaoa_bfgs_polish(const qaoa_problem_t *problem, int p, double *x,
                               int max_iter) {
  int n = 2 * p;
  size_t nn = (size_t)n;
  double *buf = malloc((6 * nn + nn * nn) * sizeof *buf);
  if (!buf) {
    return qaoa_objective(problem, p, x);
  }

  double *g = buf;
  double *g_new = g + nn;
  double *d = g_new + nn;
  double *x_new = d + nn;
  double *s = x_new + nn;
  double *y = s + nn;
  double *H = y + nn; // inverse-Hessian approximation, n x n row-major

  for (size_t i = 0; i < nn; i++) {
    for (size_t j = 0; j < nn; j++) {
      H[i * nn + j] = (i == j) ? 1.0 : 0.0;
    }
  }

  double f = qaoa_objective(problem, p, x);
  qaoa_gradient(problem, p, x, g);

  for (int iter = 0; iter < max_iter; iter++) {
    double gmax = 0.0;
    for (size_t i = 0; i < nn; i++) {
      gmax = fmax(gmax, fabs(g[i]));
    }

    if (gmax < 1e-9) {
      break;
    }

    double slope = 0.0;
    for (size_t i = 0; i < nn; i++) {
      double di = 0.0;
      for (size_t j = 0; j < nn; j++) {
        di -= H[i * nn + j] * g[j];
      }

      d[i] = di;
      slope += di * g[i];
    }

    if (slope >= 0.0) { // not a descent direction: reset to steepest descent
      for (size_t i = 0; i < nn; i++) {
        for (size_t j = 0; j < nn; j++) {
          H[i * nn + j] = (i == j) ? 1.0 : 0.0;
        }

        d[i] = -g[i];
      }

      slope = 0.0;
      for (size_t i = 0; i < nn; i++) {
        slope += d[i] * g[i];
      }
    }

    double step = 1.0;
    double f_new = f;
    int accepted = 0;
    for (int ls = 0; ls < 40; ls++) {
      for (size_t i = 0; i < nn; i++) {
        x_new[i] = x[i] + step * d[i];
      }

      f_new = qaoa_objective(problem, p, x_new);
      if (f_new <= f + 1e-4 * step * slope) {
        accepted = 1;
        break;
      }

      step *= 0.5;
    }

    if (!accepted) {
      break;
    }

    qaoa_gradient(problem, p, x_new, g_new);
    double ys = 0.0;
    for (size_t i = 0; i < nn; i++) {
      s[i] = x_new[i] - x[i];
      y[i] = g_new[i] - g[i];
      ys += y[i] * s[i];
    }

    if (ys > 1e-12) {
      // H <- (I - \rho s y^T) H (I - \rho y s^T) + \rho s s^T
      double rho = 1.0 / ys;
      double *Hy = d; // d is no longer needed this iteration
      for (size_t i = 0; i < nn; i++) {
        double acc = 0.0;
        for (size_t j = 0; j < nn; j++) {
          acc += H[i * nn + j] * y[j];
        }

        Hy[i] = acc;
      }

      double yHy = 0.0;
      for (size_t i = 0; i < nn; i++) {
        yHy += y[i] * Hy[i];
      }

      for (size_t i = 0; i < nn; i++) {
        for (size_t j = 0; j < nn; j++) {
          H[i * nn + j] += (1.0 + rho * yHy) * rho * s[i] * s[j] -
                           rho * (Hy[i] * s[j] + s[i] * Hy[j]);
        }
      }
    }

    double df = f - f_new;
    memcpy(x, x_new, nn * sizeof *x);
    memcpy(g, g_new, nn * sizeof *g);
    f = f_new;
    if (df < 1e-15 * fmax(1.0, fabs(f))) {
      break;
    }
  }

  free(buf);

  return f;
}

qaoa_result_t qaoa_run(const qaoa_problem_t *problem, int p, int n_restarts,
                       int n_sweeps, double window, uint64_t seed) {
  qaoa_result_t result = {0};

  if (!problem || p < 1 || n_restarts < 1 || n_sweeps < 1 || !(window > 0.0) ||
      !isfinite(window)) {
    return result;
  }

  int n_params = 2 * p;
  double *x = malloc((size_t)n_params * sizeof *x);
  double *best_x = malloc((size_t)n_params * sizeof *best_x);
  if (!x || !best_x) {
    free(x);
    free(best_x);

    return result;
  }

  rng_state_t rng;
  rng_seed(&rng, seed);
  double sigma = qaoa_cost_sigma(problem);
  double best_f = INFINITY;

  for (int r = 0; r < n_restarts; r++) {
    if (r == 0) {
      // Linear-ramp (annealing-inspired) start. Minimizing C is maximizing
      // -C, which flips sign of gamma.
      double sign = problem->maximize ? 1.0 : -1.0;
      for (int l = 0; l < p; l++) {
        double t = ((double)l + 0.5) / (double)p;
        x[l] = sign * (1.5 / sigma) * t;
        x[p + l] = 0.8 * (1.0 - t);
      }
    } else {
      for (int l = 0; l < p; l++) {
        x[l] = rng_uniform_range(&rng, -M_PI, M_PI);
        x[p + l] = rng_uniform_range(&rng, -0.5 * M_PI, 0.5 * M_PI);
      }
    }

    qaoa_coord_closure_t cl = {problem, p, x, 0};
    for (int sweep = 0; sweep < n_sweeps; sweep++) {
      for (int i = 0; i < n_params; i++) {
        cl.idx = i;

        x[i] = golden_section_minimize(x[i] - window, x[i] + window,
                                       qaoa_coord_eval, &cl, 1e-8);
      }
    }

    double f = qaoa_bfgs_polish(problem, p, x, 200);
    if (f < best_f) {
      best_f = f;

      memcpy(best_x, x, (size_t)n_params * sizeof *x);
    }
  }

  free(x);

  result.gamma = malloc((size_t)p * sizeof *result.gamma);
  result.beta = malloc((size_t)p * sizeof *result.beta);
  if (!result.gamma || !result.beta) {
    free(result.gamma);
    free(result.beta);
    free(best_x);

    qaoa_result_t empty = {0};

    return empty;
  }

  memcpy(result.gamma, best_x, (size_t)p * sizeof *result.gamma);
  memcpy(result.beta, best_x + p, (size_t)p * sizeof *result.beta);
  free(best_x);

  result.depth = p;

  cvector_t *psi = qaoa_prepare_state(problem, p, result.gamma, result.beta);
  if (!psi) {
    qaoa_result_free(&result);
    qaoa_result_t empty = {0};

    return empty;
  }

  result.expectation = qaoa_expect_state(problem, psi);
  result.approx_ratio = qaoa_approximation_ratio(problem, result.expectation);

  double best_cost = qaoa_problem_best_cost(problem, NULL);
  double tol = 1e-9 * fmax(1.0, fabs(best_cost));
  double p_opt = 0.0;
  double p_max = -1.0;
  for (int k = 0; k < problem->dim; k++) {
    double pk = c_abs2(psi->data[k]);
    if (fabs(problem->cost[k] - best_cost) <= tol) {
      p_opt += pk;
    }

    if (pk > p_max) {
      p_max = pk;
      result.best_bitstring = k;
    }
  }

  result.optimum_prob = p_opt;

  cvector_free(psi);

  return result;
}

void qaoa_result_free(qaoa_result_t *result) {
  if (!result) {
    return;
  }

  free(result->gamma);
  free(result->beta);

  result->gamma = NULL;
  result->beta = NULL;
}
