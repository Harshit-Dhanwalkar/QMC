/*
 * Test: Zeeman effect of a p electron from weak field to Paschen-Back
 *
 * 1. Zero field: fine-structure levels \xi/2 (x4) and -\xi (x2)
 *     D2 : D1 intensity 2 : 1
 * 2. Weak field: slopes g_J m_j with Lande factors 4/3 and 2/3
 *    (cross-checked against zeeman.h)
 * 3. Strong field: Paschen-Back energies b (m_l + 2 m_s) + \xi m_l m_s, lines
 *    gathering into normal Zeeman triplet
 * 4. Closed-form blocks agree with a numerical 6 x 6 diagonalisation of H
 * 5. Hellmann-Feynman: dE/db = <L_z + 2 S_z>; orthonormal states; zero trace;
 *    intensity sum rules
 * 6. Invalid input
 */

#include "../physics/zeeman.h"
#include "../physics/zeeman_full.h"
#include <math.h>
#include <stdio.h>

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

// H on 6-dim |m_l, m_s> basis, index (m_l + 1) * 2 + down
static void build_h(double xi, double b, double h[6][6]) {
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      h[i][j] = 0.0;
    }
  }
  for (int ml = -1; ml <= 1; ml++) {
    for (int dn = 0; dn < 2; dn++) {
      int i = (ml + 1) * 2 + dn;
      double ms = dn ? -0.5 : 0.5;

      h[i][i] = xi * ml * ms + b * (ml + 2.0 * ms);
      // (\xi / 2)(L+ S- + L- S+): up -> down lowers ms and raises ml
      if (!dn && ml < 1) { // L+ S- |ml, up> = sqrt(2 - ml(ml+1)) |ml+1, down>
        double amp = 0.5 * xi * sqrt(2.0 - ml * (ml + 1.0));

        h[(ml + 2) * 2 + 1][i] += amp;
        h[i][(ml + 2) * 2 + 1] += amp;
      }
    }
  }
}

static void jacobi(double a[6][6], double *w) {
  for (int sweep = 0; sweep < 60; sweep++) {
    double off = 0.0;

    for (int p = 0; p < 6; p++) {
      for (int q = p + 1; q < 6; q++) {
        off += a[p][q] * a[p][q];
      }
    }

    if (off < 1e-30) {
      break;
    }

    for (int p = 0; p < 6; p++) {
      for (int q = p + 1; q < 6; q++) {
        if (fabs(a[p][q]) < 1e-300) {
          continue;
        }

        double th = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        double t = (th >= 0 ? 1.0 : -1.0) / (fabs(th) + sqrt(th * th + 1.0));
        double c = 1.0 / sqrt(t * t + 1.0);
        double s = t * c;

        for (int k = 0; k < 6; k++) {
          double x = a[k][p], y = a[k][q];

          a[k][p] = c * x - s * y;
          a[k][q] = s * x + c * y;
        }

        for (int k = 0; k < 6; k++) {
          double x = a[p][k], y = a[q][k];

          a[p][k] = c * x - s * y;
          a[q][k] = s * x + c * y;
        }
      }
    }
  }

  for (int i = 0; i < 6; i++) {
    w[i] = a[i][i];
  }

  for (int i = 1; i < 6; i++) {
    double x = w[i];
    int j = i - 1;

    while (j >= 0 && w[j] > x) {
      w[j + 1] = w[j];
      j--;
    }
    w[j + 1] = x;
  }
}

static void test_zero_field(void) {
  printf("  === Test: Zero field ===\n");

  zf_state_t s[6];
  double d2 = 0.0, d1 = 0.0;

  zf_solve(2.0, 0.0, s);
  for (int k = 0; k < 4; k++) {
    check_close(s[k].energy, 1.0, 1e-12, "2P3/2 at \\xi/2");
  }
  check_close(s[4].energy, -2.0, 1e-12, "2P1/2 at -\\xi");
  check_close(s[5].energy, -2.0, 1e-12, "2P1/2 at -\\xi");

  for (int k = 0; k < 6; k++) {
    for (int ms = 0; ms < 2; ms++) {
      for (int q = -1; q <= 1; q++) {
        double i = zf_line_intensity(&s[k], ms, q);

        if (k < 4) {
          d2 += i;
        } else {
          d1 += i;
        }
      }
    }
  }
  check_close(d2 / d1, 2.0, 1e-12, "D2 : D1 total intensity");
}

static void test_weak_field_slopes(void) {
  printf("  === Test: Weak-field slopes ===\n");

  const double xi = 10.0;
  const double h = 1e-5;
  zf_state_t sp[6];
  zf_state_t sm[6];

  zf_solve(xi, h, sp);
  zf_solve(xi, -h, sm);
  double g32 = zeeman_lande_g_factor(1, 3);
  double g12 = zeeman_lande_g_factor(1, 1);

  check_close(g32, 4.0 / 3.0, 1e-12, "g(2P3/2) from zeeman.h");
  check_close(g12, 2.0 / 3.0, 1e-12, "g(2P1/2) from zeeman.h");
  for (int k = 0; k < 6; k++) {
    double slope = (sp[k].energy - sm[k].energy) / (2 * h);
    double want = (k < 4 ? g32 : g12) * 0.5 * sp[k].two_mj;

    check_close(slope, want, 1e-6, "dE/db = g_J m_j");
  }
}

static void test_paschen_back_limit(void) {
  printf("  === Test: Paschen-Back limit ===\n");

  const double xi = 1.0;
  const double b = 1000.0;
  zf_state_t s[6];

  zf_solve(xi, b, s);
  double worst = 0.0;
  double worst_line = 0.0;

  for (int k = 0; k < 6; k++) {
    double best = 1e9;

    for (int ml = -1; ml <= 1; ml++) {
      for (int dn = 0; dn < 2; dn++) {
        double ms = dn ? -0.5 : 0.5;

        if (2 * (ml + ms) == s[k].two_mj) { // same m_j
          double e = b * (ml + 2 * ms) + xi * ml * ms;

          best = fmin(best, fabs(s[k].energy - e));
        }
      }
    }

    worst = fmax(worst, best);

    // lines from both 2S states gather at -b, 0, +b
    for (int ms = 0; ms < 2; ms++) {
      for (int q = -1; q <= 1; q++) {
        if (zf_line_intensity(&s[k], ms, q) > 1e-3) {
          double e = s[k].energy - zf_s_energy(b, ms);
          double nearest = fmin(fmin(fabs(e + b), fabs(e)), fabs(e - b));

          worst_line = fmax(worst_line, nearest);
        }
      }
    }
  }

  printf(
      "  max distance to Paschen-Back energies %.2e, lines to triplet %.2e\n",
      worst, worst_line);
  check_true(worst < 5e-3, "levels approach b(m_l + 2 m_s) + \\xi m_l m_s");
  check_true(worst_line < 2.5,
             "lines lie within \\xi of normal Zeeman triplet");
}

static void test_blocks_vs_diagonalisation(void) {
  printf("  === Test: blocks versus numerical diagonalisation ===\n");

  double worst = 0.0;

  for (int c = 0; c < 8; c++) {
    double xi = 0.7 + 0.9 * c;
    double b = (double[]){0.0, 0.2, 0.9, 2.0, 5.0, -1.3, 12.0, 30.0}[c];
    double h[6][6];
    double w[6];
    double e[6];
    zf_state_t s[6];

    zf_solve(xi, b, s);
    build_h(xi, b, h);
    jacobi(h, w);
    for (int k = 0; k < 6; k++) {
      e[k] = s[k].energy;
    }

    for (int i = 1; i < 6; i++) {
      double x = e[i];
      int j = i - 1;

      while (j >= 0 && e[j] > x) {
        e[j + 1] = e[j];
        j--;
      }

      e[j + 1] = x;
    }

    for (int k = 0; k < 6; k++) {
      worst = fmax(worst, fabs(e[k] - w[k]));
    }
  }

  printf("  max |E(block) - E(Jacobi)| = %.2e\n", worst);
  check_true(worst < 1e-10,
             "closed-form blocks match the 6 x 6 diagonalisation");
}

static void test_hellmann_feynman_and_sum_rules(void) {
  printf("  === Test: Hellmann-Feynman, orthonormality, sum rules ===\n");

  const double xi = 3.0;
  double worst_hf = 0.0;
  double worst_on = 0.0;
  double worst_tr = 0.0;
  double worst_sum = 0.0;

  for (int c = 0; c < 6; c++) {
    double b = (double[]){0.1, 0.7, 1.5, 3.0, 8.0, -2.0}[c];
    double h = 1e-6;
    zf_state_t s[6];
    zf_state_t sp[6];
    zf_state_t sm[6];

    zf_solve(xi, b, s);
    zf_solve(xi, b + h, sp);
    zf_solve(xi, b - h, sm);
    double tr = 0.0;

    for (int k = 0; k < 6; k++) {
      worst_hf = fmax(worst_hf, fabs((sp[k].energy - sm[k].energy) / (2 * h) -
                                     zf_magnetic_moment(&s[k])));
      tr += s[k].energy;

      for (int l = 0; l < 6; l++) {
        double dot = 0.0;

        for (int i = 0; i < 6; i++) {
          dot += s[k].amp[i] * s[l].amp[i];
        }

        worst_on = fmax(worst_on, fabs(dot - (k == l)));
      }
    }

    worst_tr = fmax(worst_tr, fabs(tr));
    for (int ms = 0; ms < 2; ms++) {
      for (int q = -1; q <= 1; q++) {
        double sum = 0.0;

        for (int k = 0; k < 6; k++) {
          sum += zf_line_intensity(&s[k], ms, q);
        }

        worst_sum = fmax(worst_sum, fabs(sum - 1.0));
      }
    }
  }

  printf("  Hellmann-Feynman %.2e, orthonormality %.2e, trace %.2e, "
         "intensity sum rule %.2e\n",
         worst_hf, worst_on, worst_tr, worst_sum);
  check_true(worst_hf < 1e-6 && worst_on < 1e-12 && worst_tr < 1e-12 &&
                 worst_sum < 1e-12,
             "dE/db = <L_z + 2 S_z>, orthonormal, traceless, sum rule 1");
}

static void test_invalid_input(void) {
  printf("  === Test: Invalid input ===\n");

  zf_state_t s[6];

  check_true(zf_solve(NAN, 1.0, s) == -1 && zf_solve(1.0, INFINITY, s) == -1 &&
                 zf_solve(1.0, 1.0, NULL) == -1,
             "non-finite / null");
  zf_solve(1.0, 1.0, s);
  check_true(isnan(zf_line_intensity(&s[0], 0, 2)) &&
                 isnan(zf_line_intensity(&s[0], 2, 0)) &&
                 isnan(zf_magnetic_moment(NULL)),
             "bad polarisation / spin");
}

int main(void) {
  printf(" > Zeeman effect of a p electron tests:\n");

  test_zero_field();
  test_weak_field_slopes();
  test_paschen_back_limit();
  test_blocks_vs_diagonalisation();
  test_hellmann_feynman_and_sum_rules();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_zeeman_full checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
