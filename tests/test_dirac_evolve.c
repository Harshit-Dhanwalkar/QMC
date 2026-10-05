/*
 * Test: 1D Dirac time evolution (split-operator Fourier)
 *
 * 1. Invalid-input handling and closed-form step transmission
 * 2. Free evolution: probability and energy conserved to round-off; positive-
 *    and negative-energy packets move at group velocity \pm c^2 \hbar k / E
 *    (opposite directions for same k0)
 * 3. Zitterbewegung: a spinor (1, i)/\sqrt(2) wide packet at k0 = 0 oscillates
 *    as <x>(t) = -(1 - cos(2 m c^2 t/\hbar))/2 in units \hbar/(mc), excursion
 *    being exact (to packet's small momentum spread)
 * 4. Energy conserved to O(dt^2) in a smooth potential
 * 5. Klein step: transmitted probability of a packet hitting a sharp step
 *    reproduces analytic dirac_step_transmission() in normal regime and in
 *    Klein regime (where it GROWS with step height), and vanishes where step
 *    has no propagating state
 */

#include "../physics/dirac_evolve.h"
#include "complex.h"
#include "vector.h"
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

enum { N = 2048 };
static const double DX = 0.1;

static void test_invalid(void) {
  printf("  === Invalid input and closed form ===\n");

  cvector_t *u = cvector_alloc(N);
  cvector_t *l = cvector_alloc(N);
  cvector_t *short_l = cvector_alloc(N / 2);
  cvector_t *odd_u = cvector_alloc(6);
  cvector_t *odd_l = cvector_alloc(6);
  static double v[N];

  check_true(dirac_packet_1d(u, l, DX, 0.0, 1.0, 2.0, 1, 1.0, 1.0, 1.0) == 0,
             "packet: valid call succeeds");

  check_true(dirac_evolve_1d(NULL, l, NULL, DX, 0.01, 1, 1.0, 1.0, 1.0) == -1,
             "evolve(NULL upper)");
  check_true(dirac_evolve_1d(u, short_l, NULL, DX, 0.01, 1, 1.0, 1.0, 1.0) ==
                 -1,
             "evolve: mismatched lengths");
  check_true(dirac_evolve_1d(odd_u, odd_l, NULL, DX, 0.01, 1, 1.0, 1.0, 1.0) ==
                 -1,
             "evolve: length not a power of two");
  check_true(dirac_evolve_1d(u, l, NULL, 0.0, 0.01, 1, 1.0, 1.0, 1.0) == -1,
             "evolve: dx = 0");
  check_true(dirac_evolve_1d(u, l, NULL, DX, 0.0, 1, 1.0, 1.0, 1.0) == -1,
             "evolve: dt = 0");
  check_true(dirac_evolve_1d(u, l, NULL, DX, 0.01, 0, 1.0, 1.0, 1.0) == -1,
             "evolve: zero steps");
  check_true(dirac_evolve_1d(u, l, NULL, DX, 0.01, 1, -1.0, 1.0, 1.0) == -1,
             "evolve: negative mass");
  check_true(dirac_evolve_1d(u, l, NULL, DX, 0.01, 1, 1.0, 0.0, 1.0) == -1,
             "evolve: hbar = 0");
  check_true(dirac_evolve_1d(u, l, NULL, DX, 0.01, 1, 1.0, 1.0, 0.0) == -1,
             "evolve: c = 0");

  double before = dirac_norm_1d(u, l, DX);
  v[7] = NAN;
  check_true(dirac_evolve_1d(u, l, v, DX, 0.01, 1, 1.0, 1.0, 1.0) == -1,
             "evolve: non-finite potential rejected");
  check_close(dirac_norm_1d(u, l, DX), before, 0.0,
              "evolve: rejected call leaves spinor untouched");

  check_true(dirac_packet_1d(u, l, DX, 0.0, 1.0, 0.0, 1, 1.0, 1.0, 1.0) == -1,
             "packet: sigma = 0");
  check_true(dirac_packet_1d(u, l, DX, 0.0, 1.0, 2.0, 2, 1.0, 1.0, 1.0) == -1,
             "packet: branch = 2");
  check_true(dirac_packet_1d(u, l, DX, NAN, 1.0, 2.0, 1, 1.0, 1.0, 1.0) == -1,
             "packet: NaN centre");
  check_true(isnan(dirac_norm_1d(NULL, l, DX)), "norm(NULL)");
  check_true(isnan(dirac_norm_1d(u, short_l, DX)), "norm: mismatched lengths");
  check_true(isnan(dirac_position_1d(NULL, NULL, DX)), "position(NULL)");
  check_true(isnan(dirac_energy_1d(odd_u, odd_l, NULL, DX, 1.0, 1.0, 1.0)),
             "energy: length not a power of two");

  // Closed-form step transmission
  check_close(dirac_step_transmission(2.0, 0.0, 1.0, 1.0, 1.0), 1.0, 1e-14,
              "T: no step");
  check_close(dirac_step_transmission(2.0, 4.0, 1.0, 1.0, 1.0), 0.75, 1e-14,
              "T(E=2, V0=4): kappa = 3");
  check_close(dirac_step_transmission(2.0, 2.0, 1.0, 1.0, 1.0), 0.0, 0.0,
              "T: |E - V0| < m c^2 (no propagating state)");
  check_close(dirac_step_transmission(2.0, 1e9, 1.0, 1.0, 1.0),
              4.0 * sqrt(3.0) / ((1.0 + sqrt(3.0)) * (1.0 + sqrt(3.0))), 1e-6,
              "T: V0 -> infinity limit, kappa -> sqrt((E+m)/(E-m))");
  check_true(isnan(dirac_step_transmission(1.0, 0.0, 1.0, 1.0, 1.0)),
             "T: E = m c^2 rejected");
  check_true(isnan(dirac_step_transmission(2.0, 0.0, -1.0, 1.0, 1.0)),
             "T: negative mass rejected");

  cvector_free(u);
  cvector_free(l);
  cvector_free(short_l);
  cvector_free(odd_u);
  cvector_free(odd_l);
}

static void test_free_motion(void) {
  printf("  === Free evolution: conservation and group velocity ===\n");

  const double k0 = 1.0;
  const double e0 = sqrt(k0 * k0 + 1.0);
  const double dt = 0.02;
  const int steps = 500;
  cvector_t *u = cvector_alloc(N);
  cvector_t *l = cvector_alloc(N);

  for (int branch = 1; branch >= -1; branch -= 2) {
    char label[96];

    dirac_packet_1d(u, l, DX, -20.0, k0, 8.0, branch, 1.0, 1.0, 1.0);
    check_close(dirac_norm_1d(u, l, DX), 1.0, 1e-12,
                branch > 0 ? "packet(+) normalised" : "packet(-) normalised");

    double energy0 = dirac_energy_1d(u, l, NULL, DX, 1.0, 1.0, 1.0);
    snprintf(label, sizeof label, "<H> of branch %+d packet ~ %+d E_k0", branch,
             branch);
    check_close(energy0, (double)branch * e0, 5e-3, label);

    double x0 = dirac_position_1d(u, l, DX);

    if (dirac_evolve_1d(u, l, NULL, DX, dt, steps, 1.0, 1.0, 1.0) != 0) {
      check_true(0, "free evolution ran");
      continue;
    }
    double x1 = dirac_position_1d(u, l, DX);
    double v = (x1 - x0) / (dt * steps);

    snprintf(label, sizeof label, "group velocity, branch %+d", branch);
    check_close(v, (double)branch * k0 / e0, 3e-3, label);
    check_close(dirac_norm_1d(u, l, DX), 1.0, 1e-10, "free: norm conserved");
    check_close(dirac_energy_1d(u, l, NULL, DX, 1.0, 1.0, 1.0), energy0, 1e-10,
                "free: energy conserved");
  }

  cvector_free(u);
  cvector_free(l);
}

static void test_zitterbewegung(void) {
  printf("  === Zitterbewegung ===\n");

  const double dt = M_PI / 400.0; // 100 steps = \pi/4, 200 = \pi/2, 400 = \pi
  cvector_t *u = cvector_alloc(N);
  cvector_t *l = cvector_alloc(N);

  // Gaussian * (1, i)/\sqrt(2): equal-weight positive/negative energy mixture
  // whose spin precesses about z at 2 m c^2/\hbar
  dirac_packet_1d(u, l, DX, 0.0, 0.0, 5.0, 0, 1.0, 1.0, 1.0);
  for (int i = 0; i < N; i++) {
    complex_t g = u->data[i];

    u->data[i] = c_scale(g, 1.0 / sqrt(2.0));
    l->data[i] = c_scale(c_new(-g.im, g.re), 1.0 / sqrt(2.0));
  }
  check_close(dirac_norm_1d(u, l, DX), 1.0, 1e-12, "spinor (1, i) normalised");

  dirac_evolve_1d(u, l, NULL, DX, dt, 100, 1.0, 1.0, 1.0);
  check_close(dirac_position_1d(u, l, DX), -0.5, 0.03,
              "<x>(pi/4) = -(1 - cos(pi/2))/2");
  dirac_evolve_1d(u, l, NULL, DX, dt, 100, 1.0, 1.0, 1.0);
  check_close(dirac_position_1d(u, l, DX), -1.0, 0.03,
              "<x>(pi/2) = -1  (excursion hbar/(m c))");
  dirac_evolve_1d(u, l, NULL, DX, dt, 200, 1.0, 1.0, 1.0);
  check_close(dirac_position_1d(u, l, DX), 0.0, 0.01,
              "<x>(pi) = 0  (one full cycle of 2 m c^2 t/hbar)");

  cvector_free(u);
  cvector_free(l);
}

static void test_energy_in_potential(void) {
  printf("  === Energy conservation in a smooth potential ===\n");

  cvector_t *u = cvector_alloc(N);
  cvector_t *l = cvector_alloc(N);
  static double v[N];

  for (int i = 0; i < N; i++) {
    double x = ((double)i - 0.5 * N) * DX;

    v[i] = 0.8 * exp(-x * x / (2.0 * 3.0 * 3.0));
  }

  dirac_packet_1d(u, l, DX, -15.0, 1.5, 4.0, 1, 1.0, 1.0, 1.0);
  double e0 = dirac_energy_1d(u, l, v, DX, 1.0, 1.0, 1.0);

  dirac_evolve_1d(u, l, v, DX, 0.02, 1000, 1.0, 1.0, 1.0);
  check_close(dirac_norm_1d(u, l, DX), 1.0, 1e-10, "norm conserved");
  check_close(dirac_energy_1d(u, l, v, DX, 1.0, 1.0, 1.0) / e0, 1.0, 2e-4,
              "<H> conserved to O(dt^2)");

  cvector_free(u);
  cvector_free(l);
}

// Evolve a positive-energy packet (E ~ 2, k0 = sqrt(3)) into a sharp step at
// x = 0 and return the probability found at x > 0 after the packets separate.
static double klein_run(double v0, double *norm_out) {
  const double dt = 0.04;
  const int steps = 2500; // t = 100, packet reaches the step at t ~ 58
  cvector_t *u = cvector_alloc(N);
  cvector_t *l = cvector_alloc(N);
  static double v[N];

  for (int i = 0; i < N; i++) {
    double x = ((double)i - 0.5 * N) * DX;

    v[i] = x > 0.0 ? v0 : 0.0;
  }

  dirac_packet_1d(u, l, DX, -50.0, sqrt(3.0), 6.0, 1, 1.0, 1.0, 1.0);
  dirac_evolve_1d(u, l, v, DX, dt, steps, 1.0, 1.0, 1.0);

  double right = 0.0;
  for (int i = N / 2 + 1; i < N; i++) {
    right += c_abs2(u->data[i]) + c_abs2(l->data[i]);
  }
  if (norm_out) {
    *norm_out = dirac_norm_1d(u, l, DX);
  }

  cvector_free(u);
  cvector_free(l);

  return right * DX;
}

static void test_klein_step(void) {
  printf("  === Klein step (E = 2 m c^2, k0 = sqrt 3) ===\n");

  const double e0 = 2.0;
  const double heights[] = {0.0, 0.5, 4.0, 6.0, 10.0};
  double t_sim[5];

  for (int i = 0; i < 5; i++) {
    char label[96];
    double norm = 0.0;

    t_sim[i] = klein_run(heights[i], &norm);
    snprintf(label, sizeof label, "T(V0 = %4.1f)  %s", heights[i],
             heights[i] > e0 + 1.0 ? "[Klein regime]" : "[normal regime]");
    check_close(t_sim[i],
                dirac_step_transmission(e0, heights[i], 1.0, 1.0, 1.0), 0.01,
                label);
    check_close(norm, 1.0, 1e-9, "  norm conserved");
  }

  double t_evanescent = klein_run(2.0, NULL);
  check_close(t_evanescent, 0.0, 1e-3,
              "T(V0 = 2): |E - V0| < m c^2, total reflection");
  check_true(t_sim[4] > t_sim[2] && t_sim[2] > t_evanescent + 0.5,
             "Klein paradox: transmission returns and keeps growing with V0");
}

int main(void) {
  printf(" > Dirac time evolution tests\n");

  test_invalid();
  test_free_motion();
  test_zitterbewegung();
  test_energy_in_potential();
  test_klein_step();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll Dirac evolution checks passed\n");
  return 0;
}
