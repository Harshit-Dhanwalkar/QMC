/*
 * Test: quantum kicked rotor and classical standard map (kicked_rotor.h).
 *
 * 1. Quantum resonance hbar = 4 pi: <m^2> = (n K / hbar)^2 / 2 exactly
 * 2. Unitarity and the free rotor (K = 0)
 * 3. Dynamical localisation: <m^2> saturates for K = 5, hbar = 1 while the
 *    classical map keeps diffusing at D ~ K^2/2
 * 4. Correspondence: for small hbar the quantum <p^2> follows the classical
 *    diffusion before the break time
 * 5. Regular classical motion at small K stays bounded
 * 6. Invalid input
 */

#include "../physics/kicked_rotor.h"
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

enum { N = 4096, NP = 20000 };

static double g_theta[NP];
static double g_p[NP];

int main(void) {
  printf(" > Kicked rotor\n\n");
  cvector_t *psi = cvector_alloc(N);

  if (!psi) {
    return 1;
  }

  printf("  === 1. Quantum resonance (hbar = 4 pi) ===\n");
  const double hres = 4.0 * M_PI;

  krotor_init(psi);
  int done = 0;

  for (int n = 1; n <= 20; n++) {
    krotor_step(psi, 5.0, hres, 1);
    done = n;
    if (n == 1 || n == 4 || n == 20) {
      double x = n * 5.0 / hres;
      char label[48];

      snprintf(label, sizeof label, "<m^2> after %d kicks", n);
      check_close(krotor_momentum2(psi), x * x / 2.0, 1e-8, label);
    }
  }
  check_true(done == 20, "20 kicks applied");

  printf("\n  === 2. Unitarity and the free rotor ===\n");
  krotor_init(psi);
  krotor_step(psi, 7.0, 1.3, 400);
  check_close(krotor_norm(psi), 1.0, 1e-9, "norm after 400 kicks");
  krotor_init(psi);
  krotor_step(psi, 0.0, 1.0, 50);
  check_close(krotor_momentum2(psi), 0.0, 1e-20, "K = 0 stays in m = 0");

  printf("\n  === 3. Dynamical localisation (K = 5, hbar = 1) ===\n");
  krotor_init(psi);
  int ns[] = {200, 400, 800};
  double q[3];
  int t = 0;

  for (int i = 0; i < 3; i++) {
    krotor_step(psi, 5.0, 1.0, ns[i] - t);
    t = ns[i];
    q[i] = krotor_momentum2(psi);
    printf("  quantum <m^2> at n = %d: %.2f\n", t, q[i]);
  }
  double lo = fmin(q[0], fmin(q[1], q[2]));
  double hi = fmax(q[0], fmax(q[1], q[2]));

  check_true(hi < 2.0 * lo, "saturated: max < 2 min over n = 200..800");
  check_true(hi < 400.0, "stays near the localisation length");
  krotor_classical_init(g_theta, g_p, NP, 1.0, 7);
  krotor_classical_step(g_theta, g_p, NP, 5.0, 400);
  double cl = krotor_classical_p2(g_p, NP);

  printf("  classical <p^2> / n at n = 400: %.2f (K^2/2 = 12.5)\n", cl / 400.0);
  check_close(cl / 400.0, 12.5, 1.25, "classical D within 10% of K^2/2");
  check_true(cl > 20.0 * q[1], "classical spread >> quantum at n = 400");

  printf("\n  === 4. Correspondence at small hbar ===\n");
  const double hb = 0.05;

  krotor_init(psi);
  krotor_step(psi, 5.0, hb, 40);
  double qp2 = hb * hb * krotor_momentum2(psi);

  krotor_classical_init(g_theta, g_p, NP, hb, 11);
  krotor_classical_step(g_theta, g_p, NP, 5.0, 40);
  double cp2 = krotor_classical_p2(g_p, NP);

  printf("  <p^2> quantum %.3f, classical %.3f\n", qp2, cp2);
  check_close(qp2 / cp2, 1.0, 0.15, "quantum / classical <p^2>");

  printf("\n  === 5. Regular classical motion (K = 0.5) ===\n");
  krotor_classical_init(g_theta, g_p, NP, 1.0, 3);
  krotor_classical_step(g_theta, g_p, NP, 0.5, 2000);
  check_true(sqrt(krotor_classical_p2(g_p, NP)) < 2.0, "rms p bounded");

  printf("\n  === 6. Invalid input ===\n");
  cvector_t *bad = cvector_alloc(100);

  check_true(krotor_init(NULL) == -1, "NULL");
  check_true(bad && krotor_init(bad) == -1, "length not a power of two");
  check_true(krotor_step(psi, 1.0, 0.0, 1) == -1, "hbar = 0");
  check_true(krotor_step(psi, NAN, 1.0, 1) == -1, "NaN K");
  check_true(krotor_step(psi, 1.0, 1.0, 0) == -1, "steps = 0");
  check_true(krotor_classical_step(NULL, g_p, NP, 1.0, 1) == -1, "NULL theta");
  check_true(isnan(krotor_classical_p2(NULL, 3)), "p2 NULL");
  cvector_free(bad);
  cvector_free(psi);

  printf("\n  %s\n", failures ? "FAILED" : "all checks passed");

  return failures ? 1 : 0;
}
