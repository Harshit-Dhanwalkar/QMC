/*
 * Test: Aharonov-Bohm ring
 *
 * 1. Levels: -2 t \cos(2 \pi (n - \phi) / N), periodic in \phi, trace 0,
 *    bounded by 2 t; sorted output
 * 2. Persistent current = -dE/d\phi (finite difference) for every filling, and
 *    closed forms for odd and even ne
 * 3. Current is odd in phi and periodic; full band carries no current; parity
 *    effect (odd vs even ne have opposite sign)
 * 4. Two-terminal transmission: unitarity T + R = 1, 0 <= T <= 1, periodic in
 *    \phi, even in phi (reciprocity)
 * 5. Symmetric ring (leads on opposite sites): T = 0 at half a flux quantum for
 *    every energy
 * 6. Hand-solvable ring N = 4 compared with the Green's function from an
 *    independent 4 x 4 formula
 * 7. Invalid input
 */

#include "../physics/ab_ring.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>

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

int main(void) {
  printf(" === Test: Levels ===\n");
  {
    double e[16];

    ar_levels_sorted(16, 0.3, 1.0, e);
    double tr = 0.0;
    double worst = 0.0;

    for (int i = 0; i < 16; i++) {
      tr += e[i];
      worst = fmax(worst, fabs(e[i]) - 2.0);
      if (i && e[i] < e[i - 1]) {
        failures++;
      }
    }

    check_close(tr, 0.0, 1e-12, "trace of levels");
    check_true(worst <= 1e-12, "bounded by 2 t");
    check_close(ar_level(16, 3, 1.3, 1.0), ar_level(16, 2, 0.3, 1.0), 1e-12,
                "\\phi -> \\phi + 1 relabels n -> n - 1");
    {
      double e2[16];
      double d = 0.0;

      ar_levels_sorted(16, 1.3, 1.0, e2);
      for (int i = 0; i < 16; i++) {
        d = fmax(d, fabs(e[i] - e2[i]));
      }

      check_true(d < 1e-12, "sorted spectrum periodic in phi");
    }

    check_close(ar_level(16, 0, 0.0, 1.0), -2.0, 1e-12,
                "ground level at zero flux");
    check_close(ar_level(16, 0, 0.5, 1.0), ar_level(16, 1, 0.5, 1.0), 1e-12,
                "degenerate at half flux");
    check_close(ar_level(8, 3, 0.37, 2.0),
                -4.0 * cos(2 * M_PI * (3 - 0.37) / 8), 1e-12,
                "formula with t = 2");
  }

  printf("  === Test: Current = -dE/d\\phi and closed forms ===\n");
  {
    const int ns[] = {8, 11, 16};

    for (int c = 0; c < 3; c++) {
      int n = ns[c];
      double worst_fd = 0.0;
      double worst_closed = 0.0;

      for (int ne = 1; ne <= n; ne++) {
        for (int k = 0; k < 9; k++) {
          // generic, 0 < \phi < 1/2, off degeneracies
          double phi = 0.0211 + 0.0523 * k + (k % 2 ? 0.0037 : 0.0);
          double h = 1e-6;
          double fd = -(ar_energy(n, ne, phi + h, 1.0) -
                        ar_energy(n, ne, phi - h, 1.0)) /
                      (2 * h);

          worst_fd = fmax(worst_fd, fabs(fd - ar_current(n, ne, phi, 1.0)));
          if (ne < n) {
            double closed;

            if (ne % 2) {
              closed = -(4 * M_PI / n) * sin(2 * M_PI * phi / n) *
                       sin(ne * M_PI / n) / sin(M_PI / n);
            } else {
              closed = (4 * M_PI / n) * sin(2 * M_PI * (0.5 - phi) / n) *
                       sin(M_PI * ne / n) / sin(M_PI / n);
            }
            worst_closed =
                fmax(worst_closed, fabs(closed - ar_current(n, ne, phi, 1.0)));
          }
        }
      }
      printf("  N=%d: max |I - (-dE/d\\phi)| = %.2e, max |I - closed form| = "
             "%.2e\n",
             n, worst_fd, worst_closed);
      check_true(worst_fd < 1e-6,
                 "current matches finite-difference derivative");
      check_true(worst_closed < 1e-10, "current matches closed form");
    }
  }

  printf("  === Test: Symmetry and parity effect ===\n");
  {
    const int n = 20;

    check_close(ar_current(n, 5, 0.23, 1.0), -ar_current(n, 5, -0.23, 1.0),
                1e-10, "odd in \\phi (odd ne)");
    check_close(ar_current(n, 6, 0.23, 1.0), -ar_current(n, 6, -0.23, 1.0),
                1e-10, "odd in \\phi (even ne)");
    check_close(ar_current(n, 5, 0.23, 1.0), ar_current(n, 5, 1.23, 1.0), 1e-10,
                "periodic in \\phi");
    check_close(ar_current(n, n, 0.31, 1.0), 0.0, 1e-10,
                "full band carries no current");
    check_close(ar_current(n, 0, 0.31, 1.0), 0.0, 1e-12,
                "empty ring carries no current");
    check_true(ar_current(n, 5, 0.1, 1.0) < 0.0 &&
                   ar_current(n, 6, 0.1, 1.0) > 0.0,
               "odd ne diamagnetic, even ne paramagnetic at small \\phi");
  }

  printf("  === Test: transmission: unitarity, periodicity, reciprocity ===\n");
  {
    double worst_u = 0.0;
    double worst_p = 0.0;
    double worst_r = 0.0;
    double tmax = 0.0;
    double tmin = 1.0;
    const int ns[] = {6, 12, 17, 24};

    for (int c = 0; c < 4; c++) {
      int n = ns[c];

      for (int ie = 1; ie < 19; ie++) {
        for (int ip = 0; ip < 7; ip++) {
          double e = -1.9 + 0.2 * ie * 0.95;
          double phi = 0.07 + 0.13 * ip;
          ar_scatter_t s, sp, sm;

          ar_scatter(n, 0, n / 2, phi, e, 1.0, &s);
          ar_scatter(n, 0, n / 2, phi + 1.0, e, 1.0, &sp);
          ar_scatter(n, 0, n / 2, -phi, e, 1.0, &sm);

          worst_u = fmax(worst_u, fabs(s.transmission + s.reflection - 1.0));
          worst_p = fmax(worst_p, fabs(s.transmission - sp.transmission));
          worst_r = fmax(worst_r, fabs(s.transmission - sm.transmission));
          tmax = fmax(tmax, s.transmission);
          tmin = fmin(tmin, s.transmission);
        }
      }
    }

    printf("  max |T+R-1| = %.2e, periodicity %.2e, reciprocity %.2e, T range "
           "[%.4f, %.4f]\n",
           worst_u, worst_p, worst_r, tmin, tmax);
    check_true(worst_u < 1e-10, "T + R = 1");
    check_true(worst_p < 1e-10, "periodic in flux");
    check_true(worst_r < 1e-10, "T(\\phi) = T(-\\phi)");
    check_true(tmin >= -1e-12 && tmax <= 1.0 + 1e-10, "0 <= T <= 1");
  }

  printf("  === Test: destructive interference at half a flux quantum ===\n");
  {
    double worst = 0.0;
    double big = 0.0;

    for (int n = 6; n <= 30; n += 4) {
      for (int ie = 1; ie < 10; ie++) {
        double e = -1.8 + 0.4 * ie * 0.9;
        ar_scatter_t s, s0;

        ar_scatter(n, 0, n / 2, 0.5, e, 1.0, &s);
        ar_scatter(n, 0, n / 2, 0.0, e, 1.0, &s0);

        worst = fmax(worst, s.transmission);
        big = fmax(big, s0.transmission);
      }
    }

    printf("  max T(\\phi=1/2) = %.2e, max T(\\phi=0) = %.4f\n", worst, big);
    check_true(worst < 1e-12, "T(\\phi = 1/2) = 0 for opposite leads");
    check_true(big > 0.3, "T(\\phi = 0) is not zero");
  }

  printf("  === Test: N = 4 ring against an independent 4 x 4 inverse ===\n");
  {
    /* leads at sites 0 and 2, flux phase on bond 3 -> 0. Sites 1 and 3 are
     * equivalent intermediate states: eliminating them gives an effective 2 x 2
     * problem  Heff = -t^2 [ (1/E) (1 + 1) , (1/E)(1 + phase*) ; c.c., ... ] */
    double e = 0.7;
    double t = 1.0;
    double phi = 0.2;
    double gamma = sqrt(4 - e * e);

    // \exp(i * \theta) = \cos(\theta) + i * \sin(\theta)
    double theta = 2 * M_PI * phi;
    double c = cos(theta);
    double s = sin(theta);

    // a00 = e - \sigma - 2t^2/e
    // sigma = 0.5 * (e - i * \gamma)
    // a00 = 0.5*e - 2t^2/e + i * 0.5 * \gamma
    double a00_re = 0.5 * e - 2.0 * t * t / e;
    double a00_im = 0.5 * gamma;

    // a20 = -(t^2/e) * (1 + \exp(i \theta))
    // \exp(i \theta) = c + i s
    // a20 = -(t^2/e) * (1 + c) - i * (t^2/e) * s
    double a20_re = -(t * t / e) * (1.0 + c);
    double a20_im = -(t * t / e) * s;

    // det = a00 * a00 - a20 * conj(a20)
    // a00 * a00 = (a00_re^2 - a00_im^2) + i * (2 * a00_re * a00_im)
    double sq_re = a00_re * a00_re - a00_im * a00_im;
    double sq_im = 2.0 * a00_re * a00_im;

    // |a20|^2 = a20_re^2 + a20_im^2
    double a20_abs2 = a20_re * a20_re + a20_im * a20_im;

    // det = a00^2 - |a20|^2
    double det_re = sq_re - a20_abs2;
    double det_im = sq_im;

    // g20 = -a20 / det
    // |g20|^2 = |a20|^2 / |det|^2
    double det_abs2 = det_re * det_re + det_im * det_im;
    double g20_abs2 = a20_abs2 / det_abs2;

    ar_scatter_t s_out;
    ar_scatter(4, 0, 2, phi, e, t, &s_out);

    // T = \gamma^2 * |g20|^2
    check_close(s_out.transmission, gamma * gamma * g20_abs2, 1e-12,
                "T(N=4) from eliminated 2 x 2");
  }

  printf("  === Test: Invalid input ===\n");
  {
    ar_scatter_t s;
    double e[4];

    check_true(ar_scatter(3, 0, 1, 0.0, 0.0, 1.0, &s) == -1, "N < 4");
    check_true(ar_scatter(8, 0, 0, 0.0, 0.0, 1.0, &s) == -1, "same site");
    check_true(ar_scatter(8, 0, 9, 0.0, 0.0, 1.0, &s) == -1,
               "site out of range");
    check_true(ar_scatter(8, 0, 4, 0.0, 2.5, 1.0, &s) == -1,
               "energy outside the band");
    check_true(ar_scatter(8, 0, 4, NAN, 0.0, 1.0, &s) == -1, "NaN flux");
    check_true(isnan(ar_current(8, 9, 0.1, 1.0)) &&
                   isnan(ar_energy(8, -1, 0.1, 1.0)),
               "bad filling");
    check_true(ar_levels_sorted(0, 0.0, 1.0, e) == -1, "no sites");
  }

  printf(failures ? "\n  %d FAILED\n" : "\n  all checks passed\n", failures);

  return failures != 0;
}
