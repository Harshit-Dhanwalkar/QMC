/*
 * Test: quantum particle in a circular well
 *
 * 1. Bessel function against tabulated values and Wronskian-type identity
 *     J_{m-1} + J_{m+1} = 2m J_m / x
 * 2. Bessel zeros against tabulated values; zeros of J_m interlace those of
 *    J_{m+1}
 * 3. Mode counting: levels sorted, degeneracy pattern,
 *    Weyl law N(k) ~ k^2/4 - k/2
 * 4. Orthonormality of modes on a polar grid
 * 5. Node structure: n - 1 radial nodes, 2 m angular nodal lines
 * 6. Boundary condition \psi(1, \theta) = 0 and an eigen-equation check by
 * finite differences
 * 7. Invalid input
 */

#include "../physics/circular_billiard.h"
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

static cb_mode_t modes[400];

int main(void) {
  printf("  === Test: Bessel function ===\n");
  {
    check_close(cb_bessel_j(0, 1.0), 0.7651976865579666, 1e-13, "J0(1)");
    check_close(cb_bessel_j(2, 3.5), 0.4586291841943434, 1e-12, "J2(3.5)");
    check_close(cb_bessel_j(10, 5.0), 0.0014678026473, 1e-12, "J10(5)");
    check_close(cb_bessel_j(0, 0.0), 1.0, 0.0, "J0(0)");
    check_close(cb_bessel_j(3, 0.0), 0.0, 0.0, "J3(0)");
    double worst = 0.0;

    for (int m = 1; m < 20; m++) {
      for (double x = 0.3; x < 30.0; x += 1.7) {
        double lhs = cb_bessel_j(m - 1, x) + cb_bessel_j(m + 1, x),
               rhs = 2.0 * m * cb_bessel_j(m, x) / x;

        worst = fmax(worst, fabs(lhs - rhs));
      }
    }

    check_close(worst, 0.0, 1e-12, "recurrence J_{m-1} + J_{m+1} = 2m J_m / x");
  }

  printf("  === Test: Bessel zeros ===");
  {
    check_close(cb_bessel_zero(0, 1), 2.404825557695773, 1e-12, "j_{0,1}");
    check_close(cb_bessel_zero(0, 2), 5.520078110286311, 1e-12, "j_{0,2}");
    check_close(cb_bessel_zero(0, 3), 8.653727912911012, 1e-12, "j_{0,3}");
    check_close(cb_bessel_zero(1, 1), 3.831705970207512, 1e-12, "j_{1,1}");
    check_close(cb_bessel_zero(2, 1), 5.135622301840683, 1e-12, "j_{2,1}");
    check_close(cb_bessel_zero(3, 1), 6.380161895923984, 1e-12, "j_{3,1}");

    int interlace = 1;

    for (int m = 0; m < 8; m++) {
      for (int n = 1; n < 5; n++) {
        double a = cb_bessel_zero(m, n), b = cb_bessel_zero(m + 1, n),
               c = cb_bessel_zero(m, n + 1);

        if (!(a < b && b < c)) {
          interlace = 0;
        }
      }
    }

    check_true(interlace, "j_{m,n} < j_{m+1,n} < j_{m,n+1}");
  }

  printf("  === Test: mode list ===\n");
  {
    int count = cb_modes(30.0, modes, 400);

    check_true(count > 0 && count <= 400, "modes found");
    check_close((double)count / cb_weyl_count(30.0), 1.0, 0.02,
                "count / Weyl estimate at k = 30");

    int sorted = 1, dbl = 1;

    for (int i = 1; i < count; i++) {
      if (modes[i].k < modes[i - 1].k - 1e-12) {
        sorted = 0;
      }
    }

    for (int i = 0; i < count; i++) {
      if (modes[i].m == 0 && modes[i].sine) {
        dbl = 0;
      }
    }

    check_true(sorted, "sorted by k");
    check_true(dbl, "no sin mode for m = 0");
    check_close(modes[0].k, 2.404825557695773, 1e-12,
                "ground state k = j_{0,1}");
    check_true(modes[1].m == 1 && modes[2].m == 1 &&
                   fabs(modes[1].k - modes[2].k) < 1e-12,
               "first excited level is the doubly degenerate m = 1");
    check_true(cb_modes(30.0, modes, 10) > 10,
               "returns the full count when the buffer is small");
  }

  printf("  === Test: Orthonormality ===\n");
  {
    int count = cb_modes(14.0, modes, 400);
    int nr = 200;
    int nt = 400;
    double worst_diag = 0.0;
    double worst_off = 0.0;

    for (int a = 0; a < count; a += 3) {
      for (int b = a; b < count; b += 3) {
        double s = 0.0;

        for (int i = 0; i < nr; i++) {
          double r = (i + 0.5) / nr;

          for (int j = 0; j < nt; j++) {
            double th = 2.0 * M_PI * (j + 0.5) / nt;

            s += cb_mode_value(&modes[a], r, th) *
                 cb_mode_value(&modes[b], r, th) * r;
          }
        }

        s *= (1.0 / nr) * (2.0 * M_PI / nt);
        if (a == b) {
          worst_diag = fmax(worst_diag, fabs(s - 1.0));
        } else {
          worst_off = fmax(worst_off, fabs(s));
        }
      }
    }

    check_close(worst_diag, 0.0, 1e-4,
                "max |<a|a> - 1| over a subset of modes");
    check_close(worst_off, 0.0, 1e-4, "max |<a|b>| over a subset of modes");
  }

  printf("  === Test: Nodes and boundary ===");
  {
    int count = cb_modes(20.0, modes, 400);
    int ok_rad = 1;
    int ok_ang = 1;
    double bmax = 0.0;

    for (int i = 0; i < count; i++) {
      // radial nodes along a ray where the angular factor is large
      double th = modes[i].m == 0
                      ? 0.0
                      : (modes[i].sine ? M_PI / (2.0 * modes[i].m) : 0.0);
      int changes = 0;
      double prev = cb_mode_value(&modes[i], 1e-3, th);

      for (int s = 1; s < 4000; s++) {
        double v = cb_mode_value(&modes[i], 0.9999 * s / 4000.0 + 1e-3, th);

        if (v * prev < 0.0) {
          changes++;
        }
        if (v != 0.0) {
          prev = v;
        }
      }
      if (changes != modes[i].n - 1) {
        ok_rad = 0;
      }

      // angular nodal lines:
      // sign changes around a circle at r = 0.5 / k * j_{m,1}... pick r where
      // J_m(kr) != 0
      if (modes[i].m >= 1) {
        double rr = 0.999;
        int ch = 0;
        double p2 = cb_mode_value(&modes[i], rr, 1e-9);

        for (int s = 1; s <= 4000; s++) {
          double v =
              cb_mode_value(&modes[i], rr, 2.0 * M_PI * s / 4000.0 + 1e-9);

          if (v * p2 < 0.0) {
            ch++;
          }
          if (v != 0.0) {
            p2 = v;
          }
        }

        if (ch != 2 * modes[i].m) {
          ok_ang = 0;
        }
      }

      bmax = fmax(bmax, fabs(cb_mode_value(&modes[i], 1.0, 0.7)));
    }

    check_true(ok_rad, "n - 1 radial nodes for every mode with k <= 20");
    check_true(ok_ang, "2m angular nodal lines for every mode with m >= 1");
    check_close(bmax, 0.0, 1e-9, "psi vanishes on the boundary");
  }

  printf("  === Test: eigen-equation by finite differences ===\n");
  {
    int count = cb_modes(18.0, modes, 400);
    double worst = 0.0;
    double h = 1e-3;

    for (int i = 0; i < count; i += 7) {
      double x = 0.31;
      double y = -0.22;
      double f0;
      double fxp;
      double fxm;
      double fyp;
      double fym;
      double r;
      double th;

      r = hypot(x, y);
      th = atan2(y, x);
      f0 = cb_mode_value(&modes[i], r, th);
      r = hypot(x + h, y);
      th = atan2(y, x + h);
      fxp = cb_mode_value(&modes[i], r, th);
      r = hypot(x - h, y);
      th = atan2(y, x - h);
      fxm = cb_mode_value(&modes[i], r, th);
      r = hypot(x, y + h);
      th = atan2(y + h, x);
      fyp = cb_mode_value(&modes[i], r, th);
      r = hypot(x, y - h);
      th = atan2(y - h, x);
      fym = cb_mode_value(&modes[i], r, th);

      double lap = (fxp + fxm + fyp + fym - 4.0 * f0) / (h * h);
      double resid = lap + modes[i].k * modes[i].k * f0;

      worst = fmax(worst, fabs(resid) / (1.0 + modes[i].k * modes[i].k));
    }

    check_close(
        worst, 0.0, 1e-4,
        "max |Laplacian \\psi + k^2 \\psi| / (1 + k^2) at a test point");
  }

  printf("  === Test: Invalid input ===\n");
  {
    check_true(isnan(cb_bessel_j(-1, 1.0)), "negative order rejected");
    check_true(isnan(cb_bessel_j(2, -1.0)), "negative argument rejected");
    check_true(isnan(cb_bessel_zero(0, 0)), "n = 0 rejected");
    check_true(cb_modes(0.0, modes, 10) == -1, "k_max = 0 rejected");
    check_true(cb_modes(10.0, NULL, 10) == -1, "NULL buffer rejected");
    check_true(isnan(cb_mode_value(NULL, 0.5, 0.0)), "NULL mode rejected");
    check_true(isnan(cb_mode_value(&modes[0], 1.5, 0.0)), "r > 1 rejected");
  }

  if (failures) {
    printf("\n  %d check(s) FAILED\n", failures);

    return 1;
  }
  printf("\n  all checks passed\n");

  return 0;
}
