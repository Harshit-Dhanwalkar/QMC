/*
 * Test: Grover search with several marked items
 *
 * 1. State-vector evolution matches \sin^2((2k+1) \theta) for several (n, M)
 * 2. Norm is conserved
 * 3. Optimal iteration count reaches success >= 1 - M/N-ish; overshooting
 *    reduces the success (it is periodic with period \pi / \theta in 2k+1)
 * 4. M = N/4 succeeds with certainty after one iteration
 * 5. M = N/2 never changes (\theta = \pi/4: state is an eigenstate)
 * 6. Marked items remain equally likely to each other
 * 7. Agreement with library's single-target grover_search
 * 8. Invalid input
 */

#include "../physics/grover_multi.h"
#include "../physics/quantum_algorithms.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

static double amp[1 << 12];
static unsigned char marked[1 << 12];

static void mark(int n, int m) {
  for (int i = 0; i < (1 << n); i++) {
    marked[i] = 0;
  }

  // spread marked items out with a stride so they are not a contiguous block
  for (int j = 0; j < m; j++) {
    marked[(j * 7 + 3) % (1 << n)] = 1;
  }
}

// Integration-style validation of multi-target Grover behavior
int main(void) {

  // Simulated success probability follows \sin^2((2k+1) \theta)
  // and norm stays 1
  printf("  === Test: State vector vs \\sin^2((2k+1) \\theta) ===\n");
  {
    int cases[][2] = {{6, 1}, {8, 3}, {10, 5}, {10, 40}, {12, 17}};

    // Sweep through and slightly beyond optimum to cover periodic behavior
    for (unsigned c = 0; c < sizeof cases / sizeof cases[0]; c++) {
      int n = cases[c][0];
      int m = cases[c][1];
      int kopt = grover_multi_optimal_iterations(n, m);
      double worst = 0.0;
      double norm_err = 0.0;

      mark(n, m);
      grover_multi_init(amp, n);
      for (int k = 0; k <= 2 * kopt + 2; k++) {
        double s = 0.0;

        worst = fmax(worst, fabs(grover_multi_success(amp, n, marked) -
                                 grover_multi_success_exact(n, m, k)));

        for (int i = 0; i < (1 << n); i++) {
          s += amp[i] * amp[i];
        }

        norm_err = fmax(norm_err, fabs(s - 1.0));
        grover_multi_step(amp, n, marked);
      }

      char label[96];

      snprintf(label, sizeof label,
               "n=%d M=%d: max |sim - exact| over %d steps", n, m,
               2 * kopt + 3);
      check_close(worst, 0.0, 1e-12, label);
      check_close(norm_err, 0.0, 1e-10, "  norm conserved");
    }
  }

  printf("  === Test: Optimal iteration count and overshoot ===\n");
  {
    int n = 10;
    int m = 1;
    int kopt = grover_multi_optimal_iterations(n, m);

    check_true(kopt == 25, "N=1024, M=1: k* = 25");
    check_true(grover_multi_success_exact(n, m, kopt) > 0.999,
               "success > 0.999 at k*");
    check_true(grover_multi_success_exact(n, m, 2 * kopt) < 0.05,
               "doubling the iterations destroys it");
    check_true(grover_multi_success_exact(n, m, kopt) >
                   grover_multi_success_exact(n, m, kopt - 3),
               "undershoot is worse");

    check_close(grover_multi_theta(n, m), asin(1.0 / 32.0), 1e-15, "theta");
  }

  printf("  === Test: Special cases ===\n");
  {
    mark(8, 64);
    grover_multi_init(amp, 8);
    grover_multi_step(amp, 8, marked);
    check_close(grover_multi_success(amp, 8, marked), 1.0, 1e-12,
                "M = N/4: certain after one iteration");
    mark(8, 128);

    for (int i = 0; i < 128; i++) {
      marked[i] = (unsigned char)(i % 2 == 0 ? 1 : 0);
    }

    for (int i = 128; i < 256; i++) {
      marked[i] = (unsigned char)(i % 2 == 0 ? 1 : 0);
    }

    grover_multi_init(amp, 8);
    for (int k = 0; k < 5; k++) {
      grover_multi_step(amp, 8, marked);
    }

    check_close(grover_multi_success(amp, 8, marked), 0.5, 1e-12,
                "M = N/2: stuck at 1/2");
  }

  printf("  === Test: Symmetry among marked items ===\n");
  {
    int n = 8;
    int m = 5;
    double lo = 1e9;
    double hi = -1e9;

    mark(n, m);
    grover_multi_init(amp, n);
    for (int k = 0; k < 4; k++) {
      grover_multi_step(amp, n, marked);
    }

    for (int i = 0; i < (1 << n); i++) {
      if (marked[i]) {
        lo = fmin(lo, amp[i]);
        hi = fmax(hi, amp[i]);
      }
    }

    check_close(hi - lo, 0.0, 1e-13, "all marked amplitudes equal");
  }

  printf("  === Test: Agreement with grover_search ===\n");
  {
    grover_result_t g = grover_search(6, 21, 4);

    if (g.probabilities == NULL) {
      check_true(0, "grover_search returned a result");
    } else {
      check_close(g.p_target, grover_multi_success_exact(6, 1, 4), 1e-9,
                  "single-target probability");

      free(g.probabilities);
    }
  }

  printf("  === Test: Invalid input ===\n");
  {
    check_true(grover_multi_init(NULL, 4) == -1, "NULL state rejected");
    check_true(grover_multi_init(amp, 0) == -1, "zero qubits rejected");
    check_true(grover_multi_init(amp, 13) == -1, "too many qubits rejected");
    check_true(grover_multi_step(amp, 4, NULL) == -1, "NULL oracle rejected");
    check_true(isnan(grover_multi_theta(4, 0)), "M = 0 rejected");
    check_true(isnan(grover_multi_theta(4, 17)), "M > N rejected");
    check_true(grover_multi_optimal_iterations(4, 0) == -1,
               "optimal k for M = 0 rejected");
    check_true(isnan(grover_multi_success_exact(4, 1, -1)),
               "negative k rejected");
  }

  if (failures) {
    printf("\n  %d check(s) FAILED\n", failures);

    return 1;
  }
  printf("\n  all checks passed\n");

  return 0;
}
