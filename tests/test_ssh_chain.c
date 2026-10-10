/*
 * Test: Su-Schrieffer-Heeger chain
 *
 * 1. Topology: winding number 1 / 0 and Zak phase \pi / 0 on two sides of v =
 *    w
 * 2. Trivial phase: no states inside bulk gap, spectrum inside bands
 * 3. Topological phase: two edge states at E = +/- (w^2 - v^2)/w (v/w)^n,
 *    exponentially close to zero
 * 4. Edge profile: probability on A sites of cell c is (1 - r) r^c with
 *    r = (v/w)^2 = \exp(-2/xi)
 * 5. Chiral symmetry: spectrum symmetric under E -> -E; random hopping does
 *    not move edge states, a staggered potential shifts them to +/- m
 * 6. Invalid input
 */

#include "../physics/ssh_chain.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8g expected=%.8g err=%.2e\n", label, got, expected, err);
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

enum { MAXN = 40 };

static double g_bond[2 * MAXN];
static double g_eps[2 * MAXN];
static double g_e[2 * MAXN];
static double g_vec[4 * MAXN * MAXN];

int main(void) {
  printf(" > SSH chain\n\n");

  printf("  === Topology ===\n");
  check_true(ssh_winding_number(1.0, 0.5) == 0, "winding(v=1, w=0.5) = 0");
  check_true(ssh_winding_number(1.0, 1.5) == 1, "winding(v=1, w=1.5) = 1");
  check_true(ssh_winding_number(0.3, 2.0) == 1, "winding(v=0.3, w=2) = 1");
  check_true(ssh_winding_number(1.0, 1.0) == -1, "gap closing flagged");
  check_close(ssh_zak_phase(1.0, 0.5), 0.0, 1e-9, "Zak phase trivial");
  check_close(ssh_zak_phase(1.0, 1.5), M_PI, 1e-9, "Zak phase topological");
  check_close(ssh_bulk_energy(1.0, 2.0, 0.0), 3.0, 1e-12, "E(k=0) = v + w");
  check_close(ssh_bulk_energy(1.0, 2.0, M_PI), 1.0, 1e-12, "E(k=pi) = |v - w|");

  printf("  === Trivial phase (v = 1, w = 0.5, 30 cells) ===\n");
  ssh_chain_bonds(g_bond, 30, 1.0, 0.5, 0.0, 1);
  check_true(ssh_chain_solve(30, g_bond, NULL, g_e, g_vec) == 0, "solve");
  double lowest = 1e9;
  double highest = 0.0;

  for (int i = 0; i < 60; i++) {
    lowest = fmin(lowest, fabs(g_e[i]));
    highest = fmax(highest, fabs(g_e[i]));
  }

  check_true(lowest >= 0.5 - 1e-9, "no state inside gap |E| < |v - w|");
  check_true(highest <= 1.5 + 1e-9, "all states below v + w");

  printf("  === Edge states (v = 1, w = 2) ===\n");
  for (int n = 10; n <= 20; n += 5) {
    ssh_chain_bonds(g_bond, n, 1.0, 2.0, 0.0, 1);
    ssh_chain_solve(n, g_bond, NULL, g_e, NULL);
    double pred = (4.0 - 1.0) / 2.0 * pow(0.5, n);
    char label[48];

    snprintf(label, sizeof label, "edge energy, %d cells", n);
    check_close(g_e[n] / pred, 1.0, 1e-3, label);
    check_close(g_e[n - 1], -g_e[n], 1e-12, "  partner at -E");
  }

  printf("  === Edge profile ===\n");
  const int n = 30;

  ssh_chain_bonds(g_bond, n, 1.0, 2.0, 0.0, 1);
  ssh_chain_solve(n, g_bond, NULL, g_e, g_vec);
  double xi = ssh_edge_decay_length(1.0, 2.0);
  double r = exp(-2.0 / xi);

  check_close(r, 0.25, 1e-12, "r = \\exp(-2/xi) = (v/w)^2");
  for (int c = 0; c < 4; c++) {
    double p = 0.0;

    for (int s = n - 1; s <= n; s++) {
      p += g_vec[s * 2 * n + 2 * c] * g_vec[s * 2 * n + 2 * c];
    }

    char label[48];

    snprintf(label, sizeof label, "A-site probability, cell %d", c);
    check_close(p, (1.0 - r) * pow(r, c), 1e-6, label);
  }

  check_true(isnan(ssh_edge_decay_length(1.0, 0.5)),
             "no edge state in trivial phase");

  printf("  === Chiral symmetry ===\n");
  ssh_chain_bonds(g_bond, n, 1.0, 2.0, 0.5, 5);
  ssh_chain_solve(n, g_bond, NULL, g_e, NULL);

  double asym = 0.0;
  for (int i = 0; i < 2 * n; i++) {
    asym = fmax(asym, fabs(g_e[i] + g_e[2 * n - 1 - i]));
  }

  check_close(asym, 0.0, 1e-10, "spectrum symmetric under E -> -E");
  check_true(fabs(g_e[n]) < 1e-6 && fabs(g_e[n - 1]) < 1e-6,
             "edge states stay at zero under 50% hopping noise");

  ssh_chain_bonds(g_bond, n, 1.0, 2.0, 0.0, 1);
  ssh_chain_onsite(g_eps, n, 0.3);
  ssh_chain_solve(n, g_bond, g_eps, g_e, NULL);

  check_close(g_e[n], 0.3, 1e-6, "staggered potential: E = +m");
  check_close(g_e[n - 1], -0.3, 1e-6, "staggered potential: E = -m");

  printf("  === Invalid input ===\n");
  check_true(ssh_chain_bonds(NULL, 4, 1.0, 1.0, 0.0, 1) == -1, "NULL bond");
  check_true(ssh_chain_bonds(g_bond, 0, 1.0, 1.0, 0.0, 1) == -1, "no cells");
  check_true(ssh_chain_bonds(g_bond, 4, 1.0, 1.0, 1.0, 1) == -1, "noise = 1");
  check_true(ssh_chain_solve(4, NULL, NULL, g_e, NULL) == -1,
             "solve NULL bond");
  check_true(isnan(ssh_zak_phase(0.0, 0.0)), "Zak phase v = w = 0");
  check_true(isnan(ssh_bulk_energy(NAN, 1.0, 0.0)), "bulk energy NaN");

  printf("\n  %s\n", failures ? "FAILED" : "all checks passed");

  return failures ? 1 : 0;
}
