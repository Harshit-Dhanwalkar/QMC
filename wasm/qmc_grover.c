/*
 * qmc_grover.c - WebAssembly front end for Grover-search demo
 *
 * Physics: physics/grover_multi.c (state-vector Grover with several marked
 * items). This file owns database size (2^n, n <= 10), set of marked items,
 * current state and success-probability curve JavaScript plots
 */

#include "../physics/grover_multi.h"
#include <math.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { GV_MAX_N = 10, GV_MAX_ITEMS = 1 << GV_MAX_N, GV_CURVE_MAX = 400 };

static int g_n = 8;
static int g_iter = 0;
static int g_marks = 0;
static double g_amp[GV_MAX_ITEMS];
static unsigned char g_marked[GV_MAX_ITEMS];
static double g_curve[GV_CURVE_MAX];
static double g_exact[GV_CURVE_MAX];

QMC_EXPORT(qmc_gv_max_n)
int qmc_gv_max_n(void) { return GV_MAX_N; }

QMC_EXPORT(qmc_gv_curve_max)
int qmc_gv_curve_max(void) { return GV_CURVE_MAX; }

// Amplitudes of current state, length 2^n
QMC_EXPORT(qmc_gv_amp_buffer)
double *qmc_gv_amp_buffer(void) { return g_amp; }

// marked[i] = 1 for marked items
QMC_EXPORT(qmc_gv_marked_buffer)
unsigned char *qmc_gv_marked_buffer(void) { return g_marked; }

// Simulated and closed-form success probability for k = 0 .. count - 1
QMC_EXPORT(qmc_gv_curve_buffer)
double *qmc_gv_curve_buffer(void) { return g_curve; }

QMC_EXPORT(qmc_gv_exact_buffer)
double *qmc_gv_exact_buffer(void) { return g_exact; }

QMC_EXPORT(qmc_gv_n)
int qmc_gv_n(void) { return g_n; }

QMC_EXPORT(qmc_gv_iter)
int qmc_gv_iter(void) { return g_iter; }

QMC_EXPORT(qmc_gv_marks)
int qmc_gv_marks(void) { return g_marks; }

// Choose 2^n items, n in 2..10; clears marks and resets state
QMC_EXPORT(qmc_gv_set_n)
int qmc_gv_set_n(int n) {
  if (n < 2 || n > GV_MAX_N) {
    return -1;
  }

  g_n = n;
  for (int i = 0; i < GV_MAX_ITEMS; i++) {
    g_marked[i] = 0;
  }

  g_marks = 0;
  g_iter = 0;

  return grover_multi_init(g_amp, g_n);
}

// Mark (flag != 0) or unmark item i. Resets state. Refuses to mark N / 2 or
// more items
QMC_EXPORT(qmc_gv_mark)
int qmc_gv_mark(int i, int flag) {
  if (i < 0 || i >= (1 << g_n)) {
    return -1;
  }

  int want = flag != 0;

  if (want && !g_marked[i] && g_marks + 1 >= (1 << g_n) / 2) {
    return -1;
  }

  if (g_marked[i] != want) {
    g_marks += want ? 1 : -1;
    g_marked[i] = (unsigned char)want;
  }

  g_iter = 0;

  return grover_multi_init(g_amp, g_n);
}

QMC_EXPORT(qmc_gv_reset)
int qmc_gv_reset(void) {
  g_iter = 0;

  return grover_multi_init(g_amp, g_n);
}

// Apply one Grover iteration. Needs at least one marked item
QMC_EXPORT(qmc_gv_step)
int qmc_gv_step(void) {
  if (g_marks < 1) {
    return -1;
  }

  g_iter++;

  return grover_multi_step(g_amp, g_n, g_marked);
}

// Reset and apply k iterations
QMC_EXPORT(qmc_gv_goto)
int qmc_gv_goto(int k) {
  if (k < 0 || k > 100000 || g_marks < 1) {
    return -1;
  }

  grover_multi_init(g_amp, g_n);

  g_iter = 0;
  for (int i = 0; i < k; i++) {
    qmc_gv_step();
  }

  return 0;
}

QMC_EXPORT(qmc_gv_success)
double qmc_gv_success(void) {
  return grover_multi_success(g_amp, g_n, g_marked);
}

QMC_EXPORT(qmc_gv_theta)
double qmc_gv_theta(void) {
  return g_marks < 1 ? NAN : grover_multi_theta(g_n, g_marks);
}

QMC_EXPORT(qmc_gv_kopt)
int qmc_gv_kopt(void) {
  return g_marks < 1 ? -1 : grover_multi_optimal_iterations(g_n, g_marks);
}

// Closed-form success \sin^2((2k + 1) \theta) at current iteration
QMC_EXPORT(qmc_gv_exact)
double qmc_gv_exact(void) {
  return grover_multi_success_exact(g_n, g_marks < 1 ? 1 : g_marks, g_iter);
}

// Fill both curve buffers for k = 0 .. count - 1 (a fresh simulation; current
// state is left untouched)
QMC_EXPORT(qmc_gv_curve)
int qmc_gv_curve(int count) {
  static double scratch[GV_MAX_ITEMS];

  if (count < 1 || count > GV_CURVE_MAX || g_marks < 1) {
    return -1;
  }

  grover_multi_init(scratch, g_n);
  for (int k = 0; k < count; k++) {
    g_curve[k] = grover_multi_success(scratch, g_n, g_marked);
    g_exact[k] = grover_multi_success_exact(g_n, g_marks, k);

    grover_multi_step(scratch, g_n, g_marked);
  }

  return 0;
}
