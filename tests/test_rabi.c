/*
 * Test: two-level Rabi oscillations
 *
 * 1. rabi_evolve_exact reproduces rabi_excited_probability at every sampled
 *    time, both on resonance (Delta=0) and detuned
 * 2. Resonant pi-pulse: starting in the ground state, evolving to
 *    t=pi/Omega gives complete population inversion
 * 3. Unitarity: |\psi[0]|^2 + |\psi[1]|^2 stays 1 throughout
 * 4. Numerical cross-check via crank_nicolson_step, isolated from the exact
 *    propagator tests
 */

#include "../core/complex.h"
#include "../core/vector.h"
#include "../physics/rabi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.6f expected=%.6f err=%.2e\n", label, got, expected, err);
  if (!(err <= tol)) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void test_exact_matches_formula(void) {
  printf("  === Test: exact propagator matches closed-form probability ===\n");

  double tol = 1e-12;
  double cases[3][2] = {{1.0, 0.0}, {1.0, 0.5}, {2.0, 3.0}}; // {\Omega, \Delta}

  for (int c = 0; c < 3; c++) {
    double Omega = cases[c][0], Delta = cases[c][1];
    for (double t = 0.1; t < 5.0; t += 0.7) {
      cvector_t *psi = cvector_alloc(2);
      psi->data[0] = c_real(1.0);
      psi->data[1] = c_zero();
      rabi_evolve_exact(psi, t, Omega, Delta);

      double P_e_evolved = c_abs2(psi->data[1]);
      double P_e_formula = rabi_excited_probability(t, Omega, Delta);

      char label[64];
      snprintf(label, sizeof label, "\\Omega=%.1f \\Delta=%.1f t=%.2f P_e",
               Omega, Delta, t);
      check_close(P_e_evolved, P_e_formula, tol, label);

      cvector_free(psi);
    }
  }
}

static void test_resonant_pi_pulse(void) {
  printf("  === Test: resonant pi-pulse (full inversion) ===\n");

  double Omega = 2.0;
  double Delta = 0.0;
  double t_pi = M_PI / Omega;

  cvector_t *psi = cvector_alloc(2);
  psi->data[0] = c_real(1.0);
  psi->data[1] = c_zero();
  rabi_evolve_exact(psi, t_pi, Omega, Delta);

  double P_ground = c_abs2(psi->data[0]);
  double P_excited = c_abs2(psi->data[1]);
  cvector_free(psi);

  printf("  pi-pulse at t=%.4f: P_ground=%.6f P_excited=%.6f\n", t_pi, P_ground,
         P_excited);

  check_close(P_excited, 1.0, 1e-10, "P_excited at pi-pulse");
  check_close(P_ground, 0.0, 1e-10, "P_ground at pi-pulse");
}

static void test_unitarity(void) {
  printf("  === Test: unitarity (norm conservation) ===\n");

  double Omega = 1.3;
  double Delta = 0.7;

  for (double t = 0.0; t < 6.0; t += 0.5) {
    cvector_t *psi = cvector_alloc(2);
    psi->data[0] = c_real(1.0);
    psi->data[1] = c_zero();
    rabi_evolve_exact(psi, t, Omega, Delta);

    double total = c_abs2(psi->data[0]) + c_abs2(psi->data[1]);

    cvector_free(psi);

    char label[32];
    snprintf(label, sizeof label, "t=%.2f norm", t);

    check_close(total, 1.0, 1e-10, label);
  }
}

static void test_numerical_cross_check(void) {
  printf("  === Test: numerical cross-check via crank_nicolson_step ===\n");

  double hbar = 1.0;
  double Omega = 1.0;
  double Delta = 0.0;
  double t_final = M_PI / Omega;
  int steps = 20000;
  double dt = t_final / steps;

  cvector_t *psi = cvector_alloc(2);
  psi->data[0] = c_real(1.0);
  psi->data[1] = c_zero();

  int rc = rabi_evolve_numerical(psi, hbar, Omega, Delta, dt, steps);
  if (rc != 0) {
    printf("  FAIL: rabi_evolve_numerical returned error %d\n", rc);
    failures++;

    cvector_free(psi);

    return;
  }

  double P_excited_numerical = c_abs2(psi->data[1]);
  double P_excited_exact = rabi_excited_probability(t_final, Omega, Delta);

  cvector_free(psi);

  printf("  numerical P_excited=%.6f vs exact=%.6f\n", P_excited_numerical,
         P_excited_exact);

  check_close(P_excited_numerical, P_excited_exact, 5e-3,
              "numerical vs exact P_excited");
}

int main(void) {
  printf(" > Two-level Rabi oscillation tests:\n");

  test_exact_matches_formula();
  test_resonant_pi_pulse();
  test_unitarity();
  test_numerical_cross_check();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll Rabi checks passed\n");

  return 0;
}
