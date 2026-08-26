/*
 * Test: core/special/hermite.c and core/special/legendre.c.
 */

#include "../core/linalg/tridiag_eigh.h"
#include "../core/special/special.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check(int cond, const char *msg) {
  if (!cond) {
    printf("  FAIL: %s\n", msg);
    failures++;
  }
}

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  if (fabs(got - expected) > tol) {
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

/* Exact Hermite roots (physicists' convention H_n)
 *
 * Returns descending (index 0 = largest), so these are compared reversed.
 */
static const double exact_roots_n5[5] = {-2.02018287045609, -0.958572464613819,
                                         0.0, 0.958572464613819,
                                         2.02018287045609};
static const double exact_roots_n10[10] = {
    -3.43615911883774,  -2.53273167423279, -1.75668364929988, -1.03661082978951,
    -0.342901327223705, 0.342901327223705, 1.03661082978951,  1.75668364929988,
    2.53273167423279,   3.43615911883774};

static void test_hermite_basic(void) {
  printf("Test: hermite() basic values and recurrence consistency\n");

  check_close(hermite(0, 1.234), 1.0, 1e-12, "H_0(x)=1 for any x");
  check_close(hermite(1, 2.0), 4.0, 1e-12, "H_1(x)=2x");
  check_close(hermite(2, 1.0), 2.0, 1e-10, "H_2(1)=4*1-2=2");
  check_close(hermite(3, 0.0), 0.0, 1e-12, "H_3(0)=0 (odd function)");

  // H_{n+1} = 2x H_n - 2n H_{n-1}, spot-checked at a generic point
  double x = 0.73;
  for (int n = 1; n < 6; n++) {
    double lhs = hermite(n + 1, x);
    double rhs = 2.0 * x * hermite(n, x) - 2.0 * n * hermite(n - 1, x);

    check_close(lhs, rhs, 1e-9, "three-term recurrence holds");
  }
}

static void test_hermite_zeros_matches_exact_roots(void) {
  printf("Test: hermite_zeros (Golub-Welsch) matches exact roots for n=5,10");

  for (int k = 0; k < 5; k++) {
    double got = hermite_zeros(5, k);
    double expected = exact_roots_n5[5 - 1 - k]; // descending vs ascending

    check_close(got, expected, 1e-9, "n=5 root matches exact value");
  }

  for (int k = 0; k < 10; k++) {
    double got = hermite_zeros(10, k);
    double expected = exact_roots_n10[10 - 1 - k];

    check_close(got, expected, 1e-9, "n=10 root matches exact value");
  }
}

static void test_hermite_zeros_are_actually_roots(void) {
  printf(
      "Test: every returned zero genuinely satisfies H_n(x)=0, for n up to 30");

  for (int n = 2; n <= 30; n += 4) {
    double *zeros = malloc((size_t)n * sizeof(double));

    for (int k = 0; k < n; k++) {
      zeros[k] = hermite_zeros(n, k);
    }

    // Sort and check strict separation
    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        check(fabs(zeros[i] - zeros[j]) > 1e-6,
              "no two returned zeros are duplicates of each other");
      }

      double Hval = hermite(n, zeros[i]);
      /* H_n grows quickly with n and |x|, so compare relative to H_n's
       * derivative scale at that point */
      double deriv = hermite_deriv(n, zeros[i]);
      double scale = fabs(deriv) > 1.0 ? fabs(deriv) : 1.0;

      check(fabs(Hval) < 1e-6 * scale,
            "H_n(returned zero) is genuinely ~0, not garbage");
    }

    free(zeros);
  }
}

static void test_hermite_zeros_all_matches_per_index(void) {
  printf("Test: hermite_zeros_all's bulk result");

  int n = 12;
  double *bulk = malloc((size_t)n * sizeof(double));
  int count = hermite_zeros_all(n, bulk);
  check(count == n, "hermite_zeros_all returns n on success");

  for (int k = 0; k < n; k++) {
    check_close(bulk[k], hermite_zeros(n, k), 1e-12,
                "hermite_zeros_all[k] == hermite_zeros(n,k)");
  }

  check(hermite_zeros_all(0, bulk) == 0, "n<=0 returns 0 cleanly");
  check(hermite_zeros_all(5, NULL) == 0, "NULL output array returns 0 cleanly");

  free(bulk);
}

static void test_legendre_array_basic(void) {
  printf("Test: legendre_array evaluates the plain (m=0) Legendre polynomial "
         "at each point");

  const double x[4] = {-1.0, -0.5, 0.0, 1.0};
  double P[4];
  legendre_array(2, x, 4, P);

  for (int i = 0; i < 4; i++) {
    check_close(P[i], legendre(2, x[i]), 1e-12,
                "legendre_array[i] == legendre(l, x[i])");
  }

  // P_2(1)=1, P_2(-1)=1, P_2(0)=-1/2 (standard Legendre values)
  check_close(P[0], 1.0, 1e-12, "P_2(-1)=1");
  check_close(P[2], -0.5, 1e-12, "P_2(0)=-1/2");
  check_close(P[3], 1.0, 1e-12, "P_2(1)=1");

  double P0[1];
  legendre_array(3, x, 0, P0);
  legendre_array(3, NULL, 4, P0);
  legendre_array(3, x, 4, NULL);
}

int main(void) {
  printf("Tests (hermite.c, legendre.c)\n");

  test_hermite_basic();
  test_hermite_zeros_matches_exact_roots();
  test_hermite_zeros_are_actually_roots();
  test_hermite_zeros_all_matches_per_index();
  test_legendre_array_basic();

  if (failures == 0) {
    printf("\nAll Tests passed.\n");
    return 0;
  } else {
    printf("\n%d Test(s) FAILED.\n", failures);
    return 1;
  }
}
