/*
 * Test: wave-packet revivals in an infinite square well
 *
 * 1. Coefficients: no truncation leak for a moderate packet, <x> = x0, <E>
 *    matches analytic estimate
 * 2. Norm conserved at arbitrary times
 * 3. Full revival: \psi(x, 2 \pi) = \psi(x, 0)
 * 4. Half revival: \psi(x, \pi) = -\psi(1 - x, 0)
 * 5. Quarter revival: \psi(x, \pi/2) = (1 - i) / 2\psi_0(x)
 *                                    - (1 + i)/2 \psi_0(1 - x)
 * 6. Return probability |<\psi_0|\psi(t)>|^2 is 1 at 2 \pi, 0 at \pi for a
 *    packet away from centre, and 1/2 at \pi/2
 * 7. Mean position of mirror image is 1 - x0
 * 8. Invalid input
 */

#include "../physics/box_revival.h"
#include "complex.h"
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

enum { NM = 200, NX = 401 };

static double max_diff(const complex_t *a, const complex_t *b, int n) {
  double d = 0.0;

  for (int i = 0; i < n; i++) {
    d = fmax(d, hypot(a[i].re - b[i].re, a[i].im - b[i].im));
  }

  return d;
}

int main(void) {
  complex_t c[NM];
  complex_t psi0[NX];
  complex_t psi[NX];
  complex_t mirror[NX];
  double x[NX];
  double xm[NX];

  for (int j = 0; j < NX; j++) {
    x[j] = j / (double)(NX - 1);
    xm[j] = 1.0 - x[j];
  }

  const double x0 = 0.3;
  const double sigma = 0.04;
  const double n0 = 25.0;
  double kept = box_gaussian_coefficients(c, NM, x0, sigma, n0);

  printf("  --- coefficients ---\n");
  check_close(kept, 1.0, 1e-6, "weight inside the first 200 modes");
  box_wavefunction(c, NM, 0.0, x, NX, psi0);
  box_wavefunction(c, NM, 0.0, xm, NX, mirror);
  check_close(box_mean_position(x, psi0, NX), x0, 2e-3, "<x> = x0");
  {
    // free Gaussian: <k^2> = k0^2 + 1 / (4 \sigma^2) and E_n = n^2 = k^2 /
    // \pi^2
    double k0 = M_PI * n0,
           e_est = (k0 * k0 + 1.0 / (4.0 * sigma * sigma)) / (M_PI * M_PI);

    check_close(box_mean_energy(c, NM) / e_est, 1.0, 2e-3,
                "<E> / ((k0^2 + 1 / 4 \\sigma^2) / \\pi^2)");
  }

  printf("  --- norm ---\n");
  {
    double worst = 0.0;

    for (int i = 0; i < 12; i++) {
      double s = 0.0;
      double t = 0.37 + 0.91 * i;

      box_wavefunction(c, NM, t, x, NX, psi);
      for (int j = 0; j + 1 < NX; j++) {
        s += 0.5 * (x[j + 1] - x[j]) * (c_abs2(psi[j]) + c_abs2(psi[j + 1]));
      }

      worst = fmax(worst, fabs(s - 1.0));
    }

    check_close(worst, 0.0, 1e-4, "norm at 12 times");
  }

  printf("  --- revivals ---\n");
  {
    box_wavefunction(c, NM, BOX_REVIVAL_TIME, x, NX, psi);
    check_close(max_diff(psi, psi0, NX), 0.0, 1e-9,
                "t = 2 \\pi: \\psi = \\psi_0");

    box_wavefunction(c, NM, M_PI, x, NX, psi);
    for (int j = 0; j < NX; j++) {
      mirror[j] = c_neg(mirror[j]);
    }
    check_close(max_diff(psi, mirror, NX), 0.0, 1e-9,
                "t = \\pi: \\psi = -\\psi_0(1 - x)");

    box_wavefunction(c, NM, 0.0, xm, NX, mirror);
    box_wavefunction(c, NM, M_PI / 2.0, x, NX, psi);

    complex_t a = {0.5, -0.5};
    complex_t b = {-0.5, -0.5};
    complex_t cat[NX];

    for (int j = 0; j < NX; j++) {
      cat[j] = c_add(c_mul(a, psi0[j]), c_mul(b, mirror[j]));
    }

    check_close(max_diff(psi, cat, NX), 0.0, 1e-9,
                "t = pi/2: ((1-i) psi_0(x) - (1+i) psi_0(1-x)) / 2");
    check_close(box_mean_position(x, mirror, NX), 1.0 - x0, 2e-3,
                "<x> of the mirror image = 1 - x0");
  }

  printf("  --- return probability ---\n");
  {
    check_close(c_abs2(box_overlap(c, NM, BOX_REVIVAL_TIME)), 1.0, 1e-12,
                "t = 2 \\pi");
    check_close(c_abs2(box_overlap(c, NM, M_PI)), 0.0, 1e-4,
                "t = \\pi (packet left of centre: mirror is elsewhere)");
    check_close(c_abs2(box_overlap(c, NM, M_PI / 2.0)), 0.5, 1e-3,
                "t = \\pi/2");
    check_close(c_abs2(box_overlap(c, NM, 3.0 * M_PI / 2.0)), 0.5, 1e-3,
                "t = 3 \\pi/2");
    check_close(c_abs2(box_overlap(c, NM, 0.0)), 1.0, 1e-12, "t = 0");
  }

  printf("  --- density map equals box_wavefunction ---\n");
  {
    static double rho[9 * NX];
    double worst = 0.0;

    check_true(box_density_map(c, NM, x, NX, BOX_REVIVAL_TIME, 9, rho) == 0,
               "map computed");
    for (int i = 0; i < 9; i++) {
      box_wavefunction(c, NM, BOX_REVIVAL_TIME * i / 8.0, x, NX, psi);

      for (int j = 0; j < NX; j++) {
        worst = fmax(worst, fabs(rho[i * NX + j] - c_abs2(psi[j])));
      }
    }

    check_close(worst, 0.0, 1e-6,
                "max |map - |\\psi|^2| (modes below 1e-8 are skipped)");
    check_true(box_density_map(c, NM, x, NX, 1.0, 1, rho) == -1,
               "single time step rejected");
  }

  printf("  --- invalid input ---\n");
  {
    complex_t t[4];
    const double xx[2] = {0.2, 0.4};

    check_true(box_gaussian_coefficients(NULL, 4, 0.5, 0.05, 1.0) < 0.0,
               "NULL rejected");
    check_true(box_gaussian_coefficients(t, 4, 1.2, 0.05, 1.0) < 0.0,
               "x0 outside the well rejected");
    check_true(box_gaussian_coefficients(t, 4, 0.5, 0.0, 1.0) < 0.0,
               "zero width rejected");
    check_true(box_gaussian_coefficients(t, 9999, 0.5, 0.05, 1.0) < 0.0,
               "too many modes rejected");
    check_true(box_wavefunction(c, NM, NAN, xx, 2, t) == -1,
               "NaN time rejected");
    check_true(isnan(box_overlap(NULL, 3, 0.0).re), "NULL overlap is NaN");
    check_true(isnan(box_mean_energy(NULL, 3)), "NULL energy is NaN");
  }

  if (failures) {
    printf("\n  %d check(s) FAILED\n", failures);

    return 1;
  }
  printf("\n  all checks passed\n");

  return 0;
}
