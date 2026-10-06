/*
 * Test: single-qubit Bloch vector and driven, damped qubit (bloch.h)
 *
 * 1. Conversions: density <-> Bloch vector, invalid input
 * 2. Closed system: precession dv/dt = b x v about b = (\Omega, 0, \Delta)
 *    (Rodrigues rotation), Rabi formula for P_e, |v| conserved
 * 3. Damping: T1 decay of an excited qubit, transverse decay rate
 *    \gamma1/2 + gamma_phi, pure dephasing leaves z alone
 * 4. Driven + damped steady state: \rho_{ee} = (\Omega^2/4) /
 *    (\Delta^2 + \gamma^2/4 + \Omega^2/2)
 */

#include "../physics/bloch.h"
#include "../physics/rabi.h"
#include "matrix.h"
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

static double norm3(const double v[3]) {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static void test_conversions(void) {
  printf("  === Conversions and input checks ===\n");

  const double in[3] = {0.3, -0.4, 0.5};
  double out[3];
  cmatrix_t *rho = bloch_density_from_vector(in);

  check_true(rho != NULL, "density from vector");
  if (rho) {
    check_true(bloch_vector_from_density(rho, out) == 0, "vector from density");
    check_close(out[0], in[0], 1e-15, "x round trip");
    check_close(out[1], in[1], 1e-15, "y round trip");
    check_close(out[2], in[2], 1e-15, "z round trip");
    check_close(CMAT(rho, 0, 0).re + CMAT(rho, 1, 1).re, 1.0, 1e-15,
                "unit trace");
    cmatrix_free(rho);
  }

  const double ground[3] = {0.0, 0.0, 1.0};
  rho = bloch_density_from_vector(ground);
  if (rho) {
    check_close(CMAT(rho, 0, 0).re, 1.0, 0.0, "north pole = |g><g|");
    cmatrix_free(rho);
  }

  const double too_long[3] = {1.0, 1.0, 0.0};
  const double nan_v[3] = {NAN, 0.0, 0.0};
  check_true(bloch_density_from_vector(too_long) == NULL, "|v| > 1 rejected");
  check_true(bloch_density_from_vector(nan_v) == NULL, "NaN rejected");
  check_true(bloch_vector_from_density(NULL, out) == -1, "NULL rho");

  double v[3] = {0.0, 0.0, 1.0};
  check_true(bloch_evolve(NULL, 1, 0, 0, 0, 0.01, 1) == -1, "evolve: NULL v");
  check_true(bloch_evolve(v, NAN, 0, 0, 0, 0.01, 1) == -1,
             "evolve: NaN \\omega");
  check_true(bloch_evolve(v, 1, 0, -1, 0, 0.01, 1) == -1,
             "evolve: \\gamma1 < 0");
  check_true(bloch_evolve(v, 1, 0, 0, -1, 0.01, 1) == -1,
             "evolve: \\gamma_{\\phi} < 0");
  check_true(bloch_evolve(v, 1, 0, 0, 0, 0.0, 1) == -1, "evolve: dt = 0");
  check_true(bloch_evolve(v, 1, 0, 0, 0, 0.01, 0) == -1, "evolve: steps = 0");

  double bad[3] = {2.0, 0.0, 0.0};
  check_true(bloch_evolve(bad, 1, 0, 0, 0, 0.01, 1) == -1, "evolve: |v| > 1");
  check_close(bad[0], 2.0, 0.0, "rejected call leaves v untouched");
}

static void test_precession(void) {
  printf(
      "  === Closed system: precession about b = (\\Omega, 0, \\Delta) ===\n");

  const double omega = 1.3;
  const double delta = 0.7;
  const double t = 2.4;
  const double dt = 0.002;
  const int steps = (int)lround(t / dt);

  // Rodrigues rotation of v0 about b-hat by angle |b| t
  const double v0[3] = {0.2, 0.9, -0.3};
  double b[3] = {omega, 0.0, delta};
  double bn = norm3(b);
  const double k[3] = {b[0] / bn, b[1] / bn, b[2] / bn};
  double ang = bn * t;
  double kv = k[0] * v0[0] + k[1] * v0[1] + k[2] * v0[2];
  const double kxv[3] = {k[1] * v0[2] - k[2] * v0[1],
                         k[2] * v0[0] - k[0] * v0[2],
                         k[0] * v0[1] - k[1] * v0[0]};
  double want[3];

  for (int i = 0; i < 3; i++) {
    want[i] =
        v0[i] * cos(ang) + kxv[i] * sin(ang) + k[i] * kv * (1.0 - cos(ang));
  }

  double v[3] = {v0[0], v0[1], v0[2]};
  check_true(bloch_evolve(v, omega, delta, 0.0, 0.0, dt, steps) == 0,
             "evolve ran");
  check_close(v[0], want[0], 1e-8, "x vs Rodrigues");
  check_close(v[1], want[1], 1e-8, "y vs Rodrigues");
  check_close(v[2], want[2], 1e-8, "z vs Rodrigues");
  check_close(norm3(v), norm3(v0), 1e-9, "|v| conserved");

  // Rabi oscillation from the ground state: P_e = (1 - z)/2
  double g[3] = {0.0, 0.0, 1.0};
  bloch_evolve(g, omega, delta, 0.0, 0.0, dt, steps);
  check_close(0.5 * (1.0 - g[2]), rabi_excited_probability(t, omega, delta),
              1e-8, "P_e = rabi_excited_probability");

  // pi pulse on resonance flips the qubit
  double p[3] = {0.0, 0.0, 1.0};
  bloch_evolve(p, 2.0, 0.0, 0.0, 0.0, 0.001, (int)lround(M_PI / 2.0 / 0.001));
  check_close(p[2], -1.0, 1e-5, "\\pi pulse: z = -1");
}

static void test_damping(void) {
  printf("  === Damping and dephasing ===\n");

  const double dt = 0.002;
  const double t = 1.5;
  const int steps = (int)lround(t / dt);
  const double g1 = 0.8;
  const double gp = 0.5;

  double v[3] = {0.0, 0.0, -1.0};
  bloch_evolve(v, 0.0, 0.0, g1, 0.0, dt, steps);
  check_close(v[2], 1.0 - 2.0 * exp(-g1 * t), 1e-8,
              "T1: z(t) = 1 - 2 \\exp^{-g1 t}");
  check_close(norm3(v) <= 1.0 + 1e-12, 1.0, 0.0, "T1: |v| <= 1");

  double u[3] = {1.0, 0.0, 0.0};
  bloch_evolve(u, 0.0, 0.0, g1, 0.0, dt, steps);
  check_close(u[0], exp(-0.5 * g1 * t), 1e-8, "coherence decays at g1/2");
  check_close(u[2], 1.0 - exp(-g1 * t), 1e-8, "x-state also relaxes in z");

  double w[3] = {1.0, 0.0, 0.0};
  bloch_evolve(w, 0.0, 0.0, 0.0, gp, dt, steps);
  check_close(w[0], exp(-gp * t), 1e-8, "dephasing: x = \\exp^{-g_{\\phi} t}");
  check_close(w[2], 0.0, 1e-12, "dephasing leaves z alone");

  double s[3] = {0.0, 0.7, 0.0};
  bloch_evolve(s, 0.0, 0.0, g1, gp, dt, steps);
  check_close(s[1], 0.7 * exp(-(0.5 * g1 + gp) * t), 1e-8,
              "1/T2 = g1/2 + g_{\\phi}");
}

static void test_steady_state(void) {
  printf("  === Driven, damped steady state ===\n");

  const double gamma = 1.0;
  const double cases[][2] = {{1.0, 0.5}, {2.0, 0.0}, {0.5, -1.2}};

  for (int c = 0; c < 3; c++) {
    double omega = cases[c][0];
    double delta = cases[c][1];
    double v[3] = {0.0, 0.0, 1.0};
    char label[64];

    bloch_evolve(v, omega, delta, gamma, 0.0, 0.005, 8000); // t = 40
    double want = 0.25 * omega * omega /
                  (delta * delta + 0.25 * gamma * gamma + 0.5 * omega * omega);

    snprintf(label, sizeof label, "\\rho_{ee}, \\Omega=%.1f \\Delta=%.1f",
             omega, delta);
    check_close(0.5 * (1.0 - v[2]), want, 1e-6, label);
  }
}

int main(void) {
  printf(" > Bloch vector tests\n");

  test_conversions();
  test_precession();
  test_damping();
  test_steady_state();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll Bloch checks passed\n");
  return 0;
}
