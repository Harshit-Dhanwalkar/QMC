/*
 * Test: Quantum teleportation, superdense coding, CHSH/Bell inequality.
 *
 * 1. Teleportation: for a fixed unknown state (\alpha, \beta), all 4 possible
 *    Bell-measurement outcomes (forced deterministically via u1,u2) must
 *    produce, after correction, Bob's qubit exactly matching original state.
 * 2. Superdense coding: all 4 two-bit messages must decode exactly (protocol is
 *    100% deterministic by construction).
 * 3. CHSH: correlator must match the known closed form \cos(\theta - \phi);
 *    optimal angle choice must saturate Tsirelson bound 2* \sqrt(2), and a
 *    naive/suboptimal angle choice must NOT reach it (sanity check that the
 *    module isn't trivially always returning same S).
 */

#include "../core/complex.h"
#include "../physics/quantum_info.h"
#include <math.h>
#include <stdio.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.10f expected=%.10f err=%.2e\n", label, got, expected,
         err);

  if (err > tol) {
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

static void test_teleport_all_outcomes(void) {
  printf("test_teleport_all_outcomes:\n");

  complex_t alpha = c_real(0.6);
  complex_t beta = c_new(0.7368487952023082, 0.31153467384692046);

  // NOTE: u values that force each of 4 measurement outcomes: qstate_measure
  // _qubit picks outcome 0 if u < p0, else 1. Bell measurement here gives
  // p0=0.5 for both bits, so u=0.1 forces 0 and u=0.9 forces 1.
  double u_for_0 = 0.1, u_for_1 = 0.9;
  const double u1_vals[2] = {u_for_0, u_for_1};
  const double u2_vals[2] = {u_for_0, u_for_1};

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      qi_teleport_result_t r = qi_teleport(alpha, beta, u1_vals[i], u2_vals[j]);

      char label_m[64];
      snprintf(label_m, sizeof label_m,
               "m1=%d m2=%d: measured bits match forced outcome", r.m1, r.m2);
      check_true(r.m1 == i && r.m2 == j, label_m);

      char label_a[64], label_b[64];
      snprintf(label_a, sizeof label_a, "m1=%d m2=%d: Bob \\alpha.re", i, j);
      snprintf(label_b, sizeof label_b, "m1=%d m2=%d: Bob \\beta.re", i, j);
      check_close(r.bob_alpha.re, alpha.re, 1e-9, label_a);
      check_close(r.bob_beta.re, beta.re, 1e-9, label_b);
      check_close(r.bob_beta.im, beta.im, 1e-9, "Bob \\beta.im");
    }
  }
}

static void test_superdense_all_messages(void) {
  printf("test_superdense_all_messages:\n");

  for (int b1 = 0; b1 <= 1; b1++) {
    for (int b2 = 0; b2 <= 1; b2++) {
      qi_superdense_result_t r = qi_superdense(b1, b2);

      char label[64];
      snprintf(label, sizeof label, "sent (%d,%d) decodes exactly", b1, b2);
      check_true(r.decoded_bit1 == b1 && r.decoded_bit2 == b2, label);
    }
  }
}

static void test_chsh_correlator_closed_form(void) {
  printf("test_chsh_correlator_closed_form:\n");

  const double pairs[4][2] = {
      {0.0, 0.0}, {0.3, 0.7}, {1.2, -0.5}, {M_PI / 4, -M_PI / 4}};

  for (int i = 0; i < 4; i++) {
    double theta = pairs[i][0], phi = pairs[i][1];
    double e = qi_chsh_correlator(theta, phi);
    double expected = cos(theta - phi);

    char label[64];
    snprintf(label, sizeof label, "E(%.2f,%.2f) matches cos(\\theta - \\phi)",
             theta, phi);
    check_close(e, expected, 1e-9, label);
  }
}

static void test_chsh_bell_violation(void) {
  printf("test_chsh_bell_violation:\n");

  double a = 0.0, a_prime = M_PI / 2.0;
  double b = M_PI / 4.0, b_prime = 3.0 * M_PI / 4.0;

  double S = qi_chsh_S(a, a_prime, b, b_prime);
  double tsirelson = 2.0 * sqrt(2.0);

  printf("  S (optimal angles) = %.8f   Tsirelson bound = %.8f   classical "
         "bound = 2.0\n",
         S, tsirelson);

  check_close(S, tsirelson, 1e-9,
              "optimal angles saturate the Tsirelson bound 2 * \\sqrt(2)");
  check_true(
      S > 2.0,
      "quantum S exceeds the classical (local hidden variable) bound of 2");

  double S_generic = qi_chsh_S(0.1, 0.4, 0.9, 1.3);
  printf("  S (generic angles) = %.8f\n", S_generic);
  check_true(fabs(S_generic - tsirelson) > 0.05,
             "generic angle choice does not also hit Tsirelson bound");
}

int main(void) {
  test_teleport_all_outcomes();
  test_superdense_all_messages();
  test_chsh_correlator_closed_form();
  test_chsh_bell_violation();

  if (failures == 0) {
    printf("\nAll test_quantum_info checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
