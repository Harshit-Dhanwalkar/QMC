/*
 * Test: tunnelling in a symmetric double well (double_well.h)
 *
 * 1. Eigensolver on known problems: harmonic oscillator E_n = n + 1/2 and
 *    infinite square well n^2 pi^2 / 8 L^2 (a = 2L box) on same grid
 * 2. Eigenvectors: orthonormal, parity (-1)^k, ground state nodeless, k-th
 *    state has k nodes
 * 3. Double well: lowest doublet, splitting shrinks exponentially with
 *    barrier, close to WKB estimate
 * 4. A particle prepared as (\psi_0 - \psi_1)/\sqrt(2) oscillates between wells
 *    as \cos^2(\Delta t / 2), mean position -a \cos(\Delta t)
 * 5. A Gaussian packet in left well tunnels with period 2 \pi / \Delta (to a
 *    few percent); norm and energy are conserved
 * 6. Invalid input
 */

#include "../physics/double_well.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures = 0;

static double v[DW_MAX_GRID], psi[DW_MAX_STATES * DW_MAX_GRID],
    e[DW_MAX_STATES];

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

// WKB splitting for double well: \Delta ~ (\omega / \pi) * \exp(-S)
// S barrier \int \sqrt(2 (V - E0)) between turning points
static double wkb_splitting(double b, double a) {
  double w = sqrt(8.0 * b / (a * a));
  double e0 = 0.5 * w;
  // turning points: V(x) = e0 -> 1 - x^2/a^2 = \sqrt(e0 / b)
  double xt = a * sqrt(1.0 - sqrt(e0 / b));
  double s = 0.0;
  const int m = 4000;

  for (int i = 0; i < m; i++) {
    double x = -xt + 2 * xt * (i + 0.5) / m;
    double u = 1.0 - x * x / (a * a);
    double vv = b * u * u;

    s += sqrt(fmax(0.0, 2.0 * (vv - e0))) * 2 * xt / m;
  }

  return w / M_PI * exp(-s);
}

// Harmonic and infinite-square-well spectra on discretised grid
static void test_eigensolver_known_problems(void) {
  printf("  === Test: eigensolver on known problems ===\n");

  const int n = 1000;
  const double hw = 8.0, dx = dw_spacing(n, hw);

  for (int i = 0; i < n; i++) {
    double x = -hw + (i + 1) * dx;

    v[i] = 0.5 * x * x;
  }

  dw_eigen(v, n, dx, 6, e, psi);
  for (int k = 0; k < 6; k++) {
    check_close(e[k], k + 0.5, 1e-3, "harmonic E_n (O(dx^2) discretisation)");
  }

  for (int i = 0; i < n; i++) {
    v[i] = 0.0;
  }

  dw_eigen(v, n, dx, 4, e, psi);
  for (int k = 0; k < 4; k++) {
    double len = 2.0 * (hw + dx);
    double cont = (k + 1) * (k + 1) * M_PI * M_PI / (2.0 * len * len);
    double discrete = (1.0 - cos((k + 1) * M_PI / (n + 1))) / (dx * dx);

    check_close(e[k], discrete, 1e-9,
                "box E_n equals discrete-Laplacian eigenvalue");
    check_close(e[k], cont, 0.02 * cont,
                "box E_n within 2% of n^2 pi^2 / 2 L^2");
  }
}

// Orthonormality, parity and node structure of eigenstates
static void test_eigenvector_structure(void) {
  printf("  === Test: eigenvector structure ===\n");

  const int n = 800;
  const double hw = 4.0;
  const double dx = dw_spacing(n, hw);

  dw_potential(10.0, 1.5, hw, n, v);
  dw_eigen(v, n, dx, 6, e, psi);

  double worst = 0.0;
  for (int a = 0; a < 6; a++) {
    for (int b = 0; b < 6; b++) {
      double dot = 0.0;

      for (int i = 0; i < n; i++) {
        dot += psi[a * n + i] * psi[b * n + i] * dx;
      }

      worst = fmax(worst, fabs(dot - (a == b)));
    }
  }

  check_true(worst < 1e-9, "orthonormal");

  int ok = 1, nodes_ok = 1;
  for (int k = 0; k < 6; k++) {
    double par = 0.0;
    int nodes = 0;

    for (int i = 0; i < n / 2; i++) {
      par += psi[k * n + i] * psi[k * n + n - 1 - i];
    }

    if (par * ((k % 2) ? -1 : 1) <= 0.0) {
      ok = 0;
    }

    for (int i = 1; i < n; i++) {
      if (psi[k * n + i] * psi[k * n + i - 1] < 0.0) {
        nodes++;
      }
    }

    if (nodes != k) {
      nodes_ok = 0;
    }
  }

  check_true(ok, "definite parity (-1)^k");
  check_true(nodes_ok, "k-th state has k nodes");

  int pos = 1;
  for (int i = 0; i < n; i++) {
    if (psi[i] < 0.0) {
      pos = 0;
    }
  }
  check_true(pos, "ground state nodeless and positive");
}

// Splitting shrinks exponentially with barrier and matches WKB
static void test_splitting_vs_barrier(void) {
  printf("  === Test: splitting versus barrier ===\n");

  const int n = 1600;
  const double hw = 5.0;
  const double a = 1.5;
  const double dx = dw_spacing(n, hw);
  const double barriers[] = {6.0, 10.0, 16.0, 24.0};
  double prev = 1e9;

  for (int c = 0; c < 4; c++) {
    double b = barriers[c], d;

    dw_potential(b, a, hw, n, v);
    dw_eigen(v, n, dx, 4, e, psi);
    d = e[1] - e[0];

    double wk = wkb_splitting(b, a);
    printf("  B=%.0f: Delta=%.4e, WKB=%.4e, ratio %.3f; second doublet "
           "starts at %.3f\n",
           b, d, wk, d / wk, e[2]);

    check_true(d < prev && d > 0.0, "splitting positive and shrinking with B");
    check_true(d / wk > 0.5 && d / wk < 2.0,
               "within a factor 2 of WKB estimate");
    check_true(e[2] > e[1] + 0.3, "doublets separated");
    prev = d;
  }
}

// Coherent tunnelling of (psi0 - psi1)/sqrt 2 between two wells
static void test_two_state_tunnelling(void) {
  printf("  === Test: two-state tunnelling ===\n");

  const int n = 1200;
  const double hw = 5.0;
  const double a = 1.5, dx = dw_spacing(n, hw);
  static double re[DW_MAX_GRID], im[DW_MAX_GRID];
  double c[DW_MAX_STATES] = {0};

  dw_potential(8.0, a, hw, n, v);
  dw_eigen(v, n, dx, 6, e, psi);

  c[0] = 1.0 / sqrt(2.0); // (psi0 - psi1)/sqrt 2 sits in left well
  c[1] = -1.0 / sqrt(2.0);

  double d = e[1] - e[0], worst = 0.0;
  for (int k = 0; k <= 12; k++) {
    double t = k * M_PI / d / 12.0 * 2.0; // two half-periods

    dw_evolve(psi, e, c, n, 2, t, re, im);
    worst = fmax(worst,
                 fabs(dw_prob_left(re, im, n, dx) - pow(cos(0.5 * d * t), 2)));
  }

  printf("  max |P_left - cos^2(Delta t / 2)| = %.2e\n", worst);
  check_true(worst < 1e-3,
             "P_left(t) = cos^2(Delta t / 2) up to small overlap of "
             "psi0 psi1 with barrier");

  dw_evolve(psi, e, c, n, 2, 0.0, re, im);
  check_true(dw_prob_left(re, im, n, dx) > 0.99, "starts in left well");

  dw_evolve(psi, e, c, n, 2, M_PI / d, re, im);
  check_true(dw_prob_left(re, im, n, dx) < 0.01,
             "in right well after a half period pi / Delta");

  double x0 = dw_mean_x(re, im, n, dx, hw);

  check_true(x0 > a * 0.6 && x0 < a * 1.1, "mean position moves to +a");

  dw_evolve(psi, e, c, n, 2, 2.0 * M_PI / d, re, im);
  check_true(dw_prob_left(re, im, n, dx) > 0.99,
             "back in left well after 2 pi / Delta");
}

// Gaussian packet in left well tunnels with period 2 pi / Delta
static void test_gaussian_packet(void) {
  printf("  === Test: Gaussian packet ===\n");

  const int n = 1200;
  const double hw = 5.0;
  const double a = 1.5;
  const double dx = dw_spacing(n, hw);
  static double re[DW_MAX_GRID], im[DW_MAX_GRID];
  double c[DW_MAX_STATES], cap;
  const double b = 8.0;
  const double w = sqrt(8.0 * b / (a * a));
  const double sigma = 1.0 / sqrt(2.0 * w);

  dw_potential(b, a, hw, n, v);
  dw_eigen(v, n, dx, 10, e, psi);
  dw_project(psi, n, dx, 10, hw, -a, sigma, c, &cap);
  printf("  captured weight %.4f\n", cap);
  check_true(cap > 0.97, "ten states capture packet");

  double d = e[1] - e[0];
  double tmin = -1.0;
  double pmin = 2.0;
  for (int k = 0; k <= 400; k++) {
    double t = 1.5 * M_PI / d * k / 400.0;

    dw_evolve(psi, e, c, n, 10, t, re, im);
    double p = dw_prob_left(re, im, n, dx);

    if (t < 1.4 * M_PI / d && p < pmin) {
      pmin = p;
      tmin = t;
    }
  }

  printf("  lowest left probability %.3f at t = %.3f (pi / Delta = %.3f)\n",
         pmin, tmin, M_PI / d);
  check_true(pmin < 0.1, "packet leaves left well");

  check_close(tmin * d / M_PI, 1.0, 0.1,
              "half-period time within 10% of pi / Delta");

  dw_evolve(psi, e, c, n, 10, 7.7, re, im);
  double norm = 0.0, en = 0.0, e0 = 0.0;
  for (int k = 0; k < 10; k++) {
    en += c[k] * c[k] * e[k];
    e0 += c[k] * c[k];
  }

  for (int i = 0; i < n; i++) {
    norm += (re[i] * re[i] + im[i] * im[i]) * dx;
  }

  check_close(norm, e0, 1e-9, "norm conserved");

  printf("  <E> = %.4f\n", en / e0);
}

// Bad-argument rejection paths
static void test_invalid_input(void) {
  printf("  === Test: invalid input ===\n");

  double c[2], cap;

  check_true(dw_potential(-1.0, 1.0, 4.0, 100, v) == -1 &&
                 dw_potential(1.0, 0.0, 4.0, 100, v) == -1 &&
                 dw_potential(1.0, 1.0, 4.0, 2, v) == -1,
             "bad potential parameters");

  dw_potential(5.0, 1.0, 4.0, 100, v);
  v[3] += 1.0; // break symmetry
  check_true(dw_eigen(v, 100, 0.08, 2, e, psi) == -1,
             "asymmetric potential rejected");

  dw_potential(5.0, 1.0, 4.0, 100, v);

  check_true(dw_eigen(v, 100, 0.0, 2, e, psi) == -1 &&
                 dw_eigen(v, 100, 0.08, 0, e, psi) == -1 &&
                 dw_eigen(v, 100, 0.08, 20, e, psi) == -1,
             "bad solver arguments");
  check_true(dw_project(psi, 100, 0.08, 2, 4.0, 0.0, -1.0, c, &cap) == -1,
             "bad packet width");
}

int main(void) {
  printf(" > Symmetric double well tests:\n");

  test_eigensolver_known_problems();
  test_eigenvector_structure();
  test_splitting_vs_barrier();
  test_two_state_tunnelling();
  test_gaussian_packet();
  test_invalid_input();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll double-well checks passed\n");

  return 0;
}
