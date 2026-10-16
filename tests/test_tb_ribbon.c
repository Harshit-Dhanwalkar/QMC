/*
 * Test: ribbon (strip) spectra of Haldane model (tight_binding.h,
 * tb_ribbon_bands / tb_ribbon_profile)
 *
 * 1. Graphene zigzag ribbon: flat band of edge states at E = 0 between Dirac
 *    points
 * 2. Chern phase (|M| < 3 sqrt(3) |t2 sin phi|): in-gap states exist, sit on
 *    edges, and cross a Fermi level in gap exactly once per edge, with opposite
 *    velocities on opposite edges; sign of velocity on first edge equals Chern
 *    number (checked for both signs of \phi)
 * 3. Trivial phase (|M| above limit): no states in gap
 * 4. Edge-state profile: normalised and decaying away from edge
 * 5. Invalid input
 */

#include "../physics/tight_binding.h"
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

enum { W = 20, NK = 64, MAXS = 64 };

/* Bulk band edges (valence maximum, conduction minimum) from a k grid. */
static void bulk_gap_edges(const tb_model_t *m, double *vmax, double *cmin) {
  double e[2];

  *vmax = -1e30;
  *cmin = 1e30;
  for (int i = 0; i < 36; i++) {
    for (int j = 0; j < 36; j++) {
      tb_bands(m, 2.0 * M_PI * i / 36, 2.0 * M_PI * j / 36, e);
      *vmax = fmax(*vmax, e[0]);
      *cmin = fmin(*cmin, e[1]);
    }
  }
}

/* Net crossings of Fermi level ef by in-gap states on r = 0 edge (mean_cell <
 * -0.5) as k1 runs over zone: +1 for each upward crossing
 * A state counts as an edge state when it lies inside bulk gap */
static int edge_crossings(const tb_model_t *m, double ef, double vmax,
                          double cmin, int side, int *n_ingap) {
  double e[MAXS];
  double mc[MAXS];
  double prev = NAN;
  int net = 0;

  *n_ingap = 0;
  for (int ik = 0; ik <= NK; ik++) {
    double cur = NAN;

    tb_ribbon_bands(m, 2.0 * M_PI * ik / NK, W, e, mc);
    for (int s = 0; s < 2 * W; s++) {
      if (e[s] > vmax + 1e-6 && e[s] < cmin - 1e-6) {
        (*n_ingap)++;
        if (side * mc[s] > 0.5) {
          cur = e[s];
        }
      }
    }

    if (!isnan(prev) && !isnan(cur)) {
      if (prev < ef && cur >= ef) {
        net++;
      } else if (prev >= ef && cur < ef) {
        net--;
      }
    }

    prev = cur;
  }

  return net;
}

int main(void) {
  double e[MAXS];
  double mc[MAXS];

  /* Plan:
   * 1. Verify zigzag graphene ribbon has near-zero edge modes at k1 = \pi
   * 2. In Chern phase (for both signs of \phi), verify in-gap edge states cross
   * a Fermi level with opposite chirality on opposite edges
   * 3. In a trivial phase, verify no in-gap edge channels cross gap
   * 4. Verify edge-state profile properties (normalisation/localisation)
   * 5. Verify invalid-input handling paths
   */

  // Graphene zigzag ribbon supports two edge-localised near-zero modes around
  // k1 = \pi (flat edge band feature)
  printf("  --- graphene zigzag ribbon: flat edge band ---\n");
  {
    tb_model_t *g = tb_model_graphene(1.0);
    double emin = 1e30;

    tb_ribbon_bands(g, M_PI, 30, e, mc);
    for (int s = 0; s < 60; s++) {
      emin = fmin(emin, fabs(e[s]));
    }

    check_close(emin, 0.0, 1e-3, "smallest |E| at k1 = \\pi (W = 30)");
    check_true(fabs(mc[29]) > 0.5 && fabs(mc[30]) > 0.5,
               "two zero modes sit on the edges");

    tb_model_free(g);
  }

  // Topological-phase sweep over opposite flux signs to check that edge-mode
  // propagation direction flips with \phi, while opposite edges keep opposite
  // velocities/chiralities
  const double t2 = 0.1;
  const double phis[2] = {M_PI / 2, -M_PI / 2};

  for (int p = 0; p < 2; p++) {
    printf("  --- Haldane, \\phi = %+.2f \\pi, M = 0 (Chern phase) ---\n",
           phis[p] / M_PI);

    tb_model_t *m = tb_model_haldane(1.0, t2, phis[p], 0.0);
    int chern = tb_haldane_chern(1.0, t2, phis[p], 0.0), n_ingap;
    double vmax;
    double cmin;

    bulk_gap_edges(m, &vmax, &cmin);
    printf("  bulk gap: valence max %.4f, conduction min %.4f\n", vmax, cmin);

    check_close(cmin - vmax, tb_haldane_valley_gap(t2, phis[p], 0.0, 1), 1e-2,
                "bulk gap = valley gap");

    int left = edge_crossings(m, 0.2, vmax, cmin, -1, &n_ingap);
    int right = edge_crossings(m, 0.2, vmax, cmin, +1, &n_ingap);

    check_true(n_ingap > 0, "states inside the gap");
    check_true(abs(left) == 1 && abs(right) == 1, "one crossing per edge");
    check_true(left == -right, "opposite velocities on the two edges");
    check_true(left == chern, "velocity on the first edge = Chern number");

    tb_ribbon_bands(m, M_PI - 0.4, W, e, mc);
    {
      int count = 0;

      for (int s = 0; s < 2 * W; s++) {
        if (fabs(e[s]) < 0.4 * (cmin - vmax) && fabs(mc[s]) > 0.9) {
          count++;
        }
      }

      check_true(count == 2, "two edge-localised states near middle of gap");
    }

    tb_model_free(m);
  }

  printf("  --- Haldane, M = 1 (trivial phase) ---\n");
  {
    tb_model_t *m = tb_model_haldane(1.0, t2, M_PI / 2, 1.0);
    double vmax;
    double cmin;
    int n_ingap;

    bulk_gap_edges(m, &vmax, &cmin);
    check_true(tb_haldane_chern(1.0, t2, M_PI / 2, 1.0) == 0, "Chern number 0");

    edge_crossings(m, 0.2, vmax, cmin, -1, &n_ingap);
    check_true(n_ingap == 0, "no ribbon states in gap at any k1");

    tb_model_free(m);
  }

  printf("  --- edge-state profile ---\n");
  {
    tb_model_t *m = tb_model_haldane(1.0, t2, M_PI / 2, 0.0);
    double prob[W], sum = 0.0;

    tb_ribbon_bands(m, M_PI - 0.4, W, e, mc);
    int state = -1;

    for (int s = 0; s < 2 * W; s++) {
      if (mc[s] < -0.9 && (state < 0 || fabs(e[s]) < fabs(e[state]))) {
        state = s;
      }
    }

    check_true(state >= 0, "left-edge state found");
    check_true(tb_ribbon_profile(m, M_PI - 0.4, W, state, prob) == TB_OK,
               "profile computed");

    for (int r = 0; r < W; r++) {
      sum += prob[r];
    }

    check_close(sum, 1.0, 1e-9, "profile normalised");
    check_true(prob[0] > prob[W / 2] * 100.0 && prob[0] > prob[1],
               "decays away from edge");

    tb_model_free(m);
  }

  printf("  --- invalid input ---\n");
  {
    tb_model_t *m = tb_model_haldane(1.0, t2, M_PI / 2, 0.0);
    tb_model_t *chain = tb_model_alloc(1, 1);
    double prob[4];

    check_true(tb_ribbon_bands(NULL, 0.0, 4, e, NULL) == TB_ERR_INVALID,
               "NULL model rejected");
    check_true(tb_ribbon_bands(m, 0.0, 0, e, NULL) == TB_ERR_INVALID,
               "zero width rejected");
    check_true(tb_ribbon_bands(m, NAN, 4, e, NULL) == TB_ERR_INVALID,
               "NaN k rejected");
    check_true(tb_ribbon_bands(m, 0.0, 4, NULL, NULL) == TB_ERR_INVALID,
               "NULL output rejected");
    check_true(tb_ribbon_bands(chain, 0.0, 4, e, NULL) == TB_ERR_INVALID,
               "1D model rejected");
    check_true(tb_ribbon_profile(m, 0.0, 2, 99, prob) == TB_ERR_INVALID,
               "state out of range rejected");

    tb_model_free(m);
    tb_model_free(chain);
  }

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nTB Ribbon all checks passed\n");
  return 0;
}
