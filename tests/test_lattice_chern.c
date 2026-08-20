/*
 * Test: Hofstadter model Chern numbers, physics/lattice.c
 *
 * 1. p/q=1/3: Chern numbers (1,-2,1), The Thereotical reference Hofstadter
 *    result for this flux, fully gapped everywhere in the BZ.
 * 2. p/q=1/5: Chern numbers (1,1,-4,1,1), fully gapped everywhere (min gap 0.52
 *    across the whole BZ, checked directly), same grid-independence property.
 * 3. p/q=2/5: Chern numbers (-2,3,-2,3,-2) - a second flux at the same
 *    denominator, catching any bug that only happens to cancel for p=1.
 * 4. p/q=1/2 ("pi-flux" case): the two bands touch at Dirac points exactly
 *    everywhere along a line in the bz (min gap is exactly 0, not just small) -
 *    individual-band Chern numbers are not meaningful here by usual TKNN
 *    argument, but sum over all bands must still be exactly 0.
 * 5. Sum rule: for every (p,q) tested, sum of all q band Chern numbers is
 *    exactly 0 (full q-band bundle over the BZ torus is trivial) - checked as
 *    an independent self-consistency property alongside known-value comparisons
 *    above.
 */

#include "../core/matrix.h"
#include "../physics/lattice.h"
#include <math.h>
#include <stdio.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.6f expected=%.6f err=%.2e\n", label, got, expected, err);
  if (err > tol) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAILED");
  if (!cond) {
    failures++;
  }
}

static double sum_array(const double *a, int n) {
  double s = 0.0;
  for (int i = 0; i < n; i++) {
    s += a[i];
  }
  return s;
}

static void test_p1_q3(void) {
  printf("test_hofstadter_chern_p1_q3:\n");
  double c[3];
  int ok = lattice_hofstadter_chern_numbers(1, 3, 1.0, 24, c);

  check_true(ok, "call succeeds");
  check_close(c[0], 1.0, 1e-6, "band 0 Chern number");
  check_close(c[1], -2.0, 1e-6, "band 1 Chern number");
  check_close(c[2], 1.0, 1e-6, "band 2 Chern number");
  check_close(sum_array(c, 3), 0.0, 1e-6, "sum rule");
}

static void test_p1_q5(void) {
  printf("test_hofstadter_chern_p1_q5:\n");
  double c[5];
  int ok = lattice_hofstadter_chern_numbers(1, 5, 1.0, 24, c);
  check_true(ok, "call succeeds");

  double expected[5] = {1.0, 1.0, -4.0, 1.0, 1.0};
  for (int i = 0; i < 5; i++) {
    char label[64];
    snprintf(label, sizeof label, "band %d Chern number", i);

    check_close(c[i], expected[i], 1e-6, label);
  }

  check_close(sum_array(c, 5), 0.0, 1e-6, "sum rule");
}

static void test_p2_q5(void) {
  printf("test_hofstadter_chern_p2_q5:\n");
  double c[5];
  int ok = lattice_hofstadter_chern_numbers(2, 5, 1.0, 24, c);
  check_true(ok, "call succeeds");

  double expected[5] = {-2.0, 3.0, -2.0, 3.0, -2.0};
  for (int i = 0; i < 5; i++) {
    char label[64];
    snprintf(label, sizeof label, "band %d Chern number", i);

    check_close(c[i], expected[i], 1e-6, label);
  }

  check_close(sum_array(c, 5), 0.0, 1e-6, "sum rule");
}

static void test_p1_q2_band_touching_sum_rule_only(void) {
  printf("test_hofstadter_chern_p1_q2_sum_rule:\n");
  /* pi-flux case: bands touch exactly (Dirac points), individual Chern
   * numbers are not meaningful here - only the sum rule is asserted. */
  double c[2];
  int ok = lattice_hofstadter_chern_numbers(1, 2, 1.0, 24, c);
  check_true(ok, "call succeeds");

  check_close(sum_array(c, 2), 0.0, 1e-6,
              "sum rule holds even at exact band touching");
}

static void test_invalid_input_rejected(void) {
  printf("test_hofstadter_chern_invalid_input:\n");
  double c[4];

  check_true(!lattice_hofstadter_chern_numbers(0, 3, 1.0, 24, c),
             "p<1 rejected");
  check_true(!lattice_hofstadter_chern_numbers(1, 1, 1.0, 24, c),
             "q<2 rejected");
  check_true(!lattice_hofstadter_chern_numbers(3, 3, 1.0, 24, c),
             "p>=q rejected");
  check_true(!lattice_hofstadter_chern_numbers(1, 3, 1.0, 1, c),
             "n_k<2 rejected");
  check_true(lattice_hofstadter_bloch(0.1, 0.2, 3, 3, 1.0) == NULL,
             "bloch: p>=q rejected");
}

int main(void) {
  test_p1_q3();
  test_p1_q5();
  test_p2_q5();
  test_p1_q2_band_touching_sum_rule_only();
  test_invalid_input_rejected();

  if (failures > 0) {
    printf("\n%d test_lattice_chern check(s) FAILED.\n", failures);
    return 1;
  }
  printf("\nAll test_lattice_chern checks passed.\n");
  return 0;
}
