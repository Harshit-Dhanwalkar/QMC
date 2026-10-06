/*
 * Test: Monte Carlo sampling of hydrogen orbitals
 *
 * 1. Invalid input, determinism, closed forms of R_nl and |Y_lm|^2
 * 2. Radial normalisation of R_nl
 * 3. Sampled <r> and <r^2> against exact hydrogen moments
 *      <r>   = (3 n^2 - l(l+1)) / 2
 *      <r^2> = n^2 (5 n^2 + 1 - 3 l(l+1)) / 2
 * 4. Angular statistics: <\cos^2(\theta)> of Y_l^m, and p_x / p_y / p_z shapes
 *    of real forms
 * 5. Phase: complex orbitals carry arg(\psi) = m \phi; orbitals have a sign
 *    flip across plane of a p orbital and across node of 2s
 */

#include "../physics/orbital_sample.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double *pt(const double *p, int i) { return p + 4 * (size_t)i; }

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.6g expected=%.6g err=%.2e\n", label, got, expected, err);
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

enum { NS = 200000 };

static double *draw(int n, int l, int m, int real_form,
                    unsigned long long seed) {
  double *pts = malloc(4 * (size_t)NS * sizeof *pts);

  if (!pts) {
    return NULL;
  }
  if (hydrogen_orbital_sample(n, l, m, real_form, NS, seed, pts) != 0) {
    free(pts);

    return NULL;
  }

  return pts;
}

static void test_basics(void) {
  printf("\n  === Input checks and closed forms ===\n");

  double buf[8];
  check_true(hydrogen_orbital_sample(0, 0, 0, 0, 2, 1, buf) == -1, "n = 0");
  check_true(hydrogen_orbital_sample(ORBITAL_MAX_N + 1, 0, 0, 0, 2, 1, buf) ==
                 -1,
             "n too large");
  check_true(hydrogen_orbital_sample(2, 2, 0, 0, 2, 1, buf) == -1, "l >= n");
  check_true(hydrogen_orbital_sample(2, -1, 0, 0, 2, 1, buf) == -1, "l < 0");
  check_true(hydrogen_orbital_sample(3, 1, 2, 0, 2, 1, buf) == -1, "|m| > l");
  check_true(hydrogen_orbital_sample(1, 0, 0, 0, 0, 1, buf) == -1, "count = 0");
  check_true(hydrogen_orbital_sample(1, 0, 0, 0, 2, 1, NULL) == -1,
             "out = NULL");
  check_true(isnan(hydrogen_orbital_radial_au(2, 2, 1.0)), "radial: l >= n");
  check_true(isnan(hydrogen_orbital_radial_au(1, 0, -1.0)), "radial: r < 0");
  check_true(isnan(hydrogen_orbital_angular_density(1, 2, 0, 0.3, 0.3)),
             "angular: |m| > l");

  check_close(hydrogen_orbital_radial_au(1, 0, 1.0), 2.0 * exp(-1.0), 1e-9,
              "R_10(1) = 2/e");
  check_close(hydrogen_orbital_radial_au(2, 1, 2.0),
              2.0 / (2.0 * sqrt(6.0)) * exp(-1.0), 1e-9,
              "R_21(2) = r e^(-r/2) / (2 sqrt 6)");
  check_close(hydrogen_orbital_angular_density(0, 0, 0, 1.0, 2.0),
              1.0 / (4.0 * M_PI), 1e-12, "|Y_00|^2 = 1/4pi");
  check_close(hydrogen_orbital_angular_density(1, 0, 1, 0.0, 0.0),
              3.0 / (4.0 * M_PI), 1e-12, "|Y_10(north pole)|^2 = 3/4pi");

  double a[8];
  double b[8];
  double c[8];

  hydrogen_orbital_sample(3, 1, 0, 0, 2, 42, a);
  hydrogen_orbital_sample(3, 1, 0, 0, 2, 42, b);
  hydrogen_orbital_sample(3, 1, 0, 0, 2, 43, c);

  int same = 1;
  int diff = 0;
  for (int i = 0; i < 8; i++) {
    same &= (a[i] == b[i]);
    diff |= (a[i] != c[i]);
  }

  check_true(same, "same seed reproduces sample exactly");
  check_true(diff, "different seed changes sample");
}

static void test_radial_norm(void) {
  printf("\n  === Radial normalisation ===\n");

  const int nl[][2] = {{1, 0}, {2, 1}, {3, 2}, {4, 3}, {5, 0}, {6, 2}};

  for (int k = 0; k < 6; k++) {
    int n = nl[k][0];
    int l = nl[k][1];
    double rmax = n * (2.0 * n + 16.0);
    int steps = 20000;
    double dr = rmax / steps;
    double sum = 0.0;
    char label[64];

    for (int i = 0; i < steps; i++) {
      double r = (i + 0.5) * dr;
      double R = hydrogen_orbital_radial_au(n, l, r);

      sum += R * R * r * r * dr;
    }

    snprintf(label, sizeof label, "\\int R^2 r^2 dr, n=%d l=%d", n, l);
    check_close(sum, 1.0, 1e-6, label);
  }
}

static void test_radial_moments(void) {
  printf("\n  === Sampled radial moments vs exact ===\n");

  const int nl[][2] = {{1, 0}, {2, 0}, {2, 1}, {3, 1}, {4, 3}, {5, 2}};

  for (int k = 0; k < 6; k++) {
    int n = nl[k][0];
    int l = nl[k][1];
    double *p = draw(n, l, 0, 0, 1000 + k);
    char label[64];

    if (!p) {
      check_true(0, "sampling succeeded");
      continue;
    }

    double s1 = 0.0;
    double s2 = 0.0;
    for (int i = 0; i < NS; i++) {
      double r2 = pt(p, i)[0] * pt(p, i)[0] + pt(p, i)[1] * pt(p, i)[1] +
                  pt(p, i)[2] * pt(p, i)[2];

      s1 += sqrt(r2);
      s2 += r2;
    }

    double ex1 = (3.0 * n * n - l * (l + 1.0)) / 2.0;
    double ex2 = n * n * (5.0 * n * n + 1.0 - 3.0 * l * (l + 1.0)) / 2.0;

    snprintf(label, sizeof label, "<r>   n=%d l=%d", n, l);
    check_close(s1 / NS, ex1, 0.012 * ex1, label);
    snprintf(label, sizeof label, "<r^2> n=%d l=%d", n, l);
    check_close(s2 / NS, ex2, 0.02 * ex2, label);

    free(p);
  }
}

static void test_angular(void) {
  printf("\n  === Angular statistics ===\n");

  const int lm[][2] = {{1, 0}, {1, 1}, {2, 0}, {2, 2}, {3, 1}};

  for (int k = 0; k < 5; k++) {
    int l = lm[k][0];
    int m = lm[k][1];
    double *p = draw(l + 1, l, m, 0, 2000 + k);
    char label[64];

    if (!p) {
      check_true(0, "sampling succeeded");
      continue;
    }

    double s = 0.0;
    for (int i = 0; i < NS; i++) {
      double r2 = pt(p, i)[0] * pt(p, i)[0] + pt(p, i)[1] * pt(p, i)[1] +
                  pt(p, i)[2] * pt(p, i)[2];

      s += pt(p, i)[2] * pt(p, i)[2] / r2;
    }

    double ex = (2.0 * l * (l + 1.0) - 2.0 * m * m - 1.0) /
                ((2.0 * l - 1.0) * (2.0 * l + 3.0));

    snprintf(label, sizeof label, "<\\cos^2(\\theta)>, l=%d m=%d", l, m);
    check_close(s / NS, ex, 0.006, label);

    free(p);
  }

  // Real p orbitals: density ~ (axis)^2 / r^2, so <axis^2/r^2> = 3/5 and other
  // two are 1/5. m = 0 -> z, m = +1 -> x (\cos(\phi)), m = -1 -> y (\sin(\phi))
  const int ms[3] = {0, 1, -1};
  const int axis[3] = {2, 0, 1};
  for (int k = 0; k < 3; k++) {
    double *p = draw(2, 1, ms[k], 1, 3000 + k);
    double q[3] = {0.0, 0.0, 0.0};

    if (!p) {
      check_true(0, "sampling succeeded");
      continue;
    }
    for (int i = 0; i < NS; i++) {
      double r2 = pt(p, i)[0] * pt(p, i)[0] + pt(p, i)[1] * pt(p, i)[1] +
                  pt(p, i)[2] * pt(p, i)[2];

      for (int a = 0; a < 3; a++) {
        q[a] += pt(p, i)[a] * pt(p, i)[a] / r2;
      }
    }

    char label[64];
    snprintf(label, sizeof label, "real p (m=%+d): <%c^2/r^2> = 3/5", ms[k],
             "xyz"[axis[k]]);
    check_close(q[axis[k]] / NS, 0.6, 0.006, label);

    snprintf(label, sizeof label, "real p (m=%+d): other axes = 1/5", ms[k]);
    check_close((q[0] + q[1] + q[2] - q[axis[k]]) / NS, 0.4, 0.008, label);

    free(p);
  }
}

static void test_phase(void) {
  printf("  === Phase and sign lobes ===\n");

  // 3d, m = 2 has no radial node and Y_22 ~ \exp^{2 i \phi}: phase = 2 \phi
  double *p = draw(3, 2, 2, 0, 4000);
  if (p) {
    double worst = 0.0;

    for (int i = 0; i < NS; i += 7) {
      double phi = atan2(pt(p, i)[1], pt(p, i)[0]);
      double d = fmod(pt(p, i)[3] - 2.0 * phi, 2.0 * M_PI);

      if (d > M_PI) {
        d -= 2.0 * M_PI;
      } else if (d < -M_PI) {
        d += 2.0 * M_PI;
      }
      if (fabs(d) > worst) {
        worst = fabs(d);
      }
    }

    check_close(worst, 0.0, 1e-6, "complex 3d(m=2): max |arg psi - 2 phi|");

    free(p);
  }

  // real 2p_z: sign set by hemisphere
  p = draw(2, 1, 0, 1, 4001);
  if (p) {
    int upper_phase = -1;
    int consistent = 1;

    for (int i = 0; i < NS; i++) {
      int up = pt(p, i)[2] > 0.0;
      int ph = pt(p, i)[3] > 1.0;

      if (up) {
        if (upper_phase < 0) {
          upper_phase = ph;
        }
        consistent &= (ph == upper_phase);
      } else {
        consistent &= (ph != upper_phase);
      }
    }

    check_true(consistent, "real 2p_z: phase is + above, - below xy plane");

    free(p);
  }

  // real 2s: radial node at r = 2 a0 flips sign
  p = draw(2, 0, 0, 1, 4002);
  if (p) {
    int inner = -1;
    int consistent = 1;
    int near_node = 0;

    for (int i = 0; i < NS; i++) {
      double r = sqrt(pt(p, i)[0] * pt(p, i)[0] + pt(p, i)[1] * pt(p, i)[1] +
                      pt(p, i)[2] * pt(p, i)[2]);
      int ph = pt(p, i)[3] > 1.0;

      if (fabs(r - 2.0) < 0.02) {
        near_node++;
      }
      if (fabs(r - 2.0) < 1e-3) {
        continue; // sign undefined on node itself
      }
      if (r < 2.0) {
        if (inner < 0) {
          inner = ph;
        }
        consistent &= (ph == inner);
      } else {
        consistent &= (ph != inner);
      }
    }

    check_true(consistent, "real 2s: sign flips across r = 2 a0");
    check_true(near_node < NS / 200, "real 2s: node is nearly empty");

    free(p);
  }
}

int main(void) {
  printf(" > Hydrogen orbital sampling tests\n");

  test_basics();
  test_radial_norm();
  test_radial_moments();
  test_angular();
  test_phase();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll orbital sampling checks passed\n");
  return 0;
}
