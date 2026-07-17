/*
Test: angular momentum coupling (couple_states / couple_allowed_J).

1. Two spin-1/2's -> triplet (J=1) + singlet (J=0), checked against
  exact coefficients (in Sakurai angular-momenta table).
2. General sanity checks for less trivial case (j1=1, j2=1/2, i.e.
   l(x)s structure needed for spin-orbit coupling): normalization of
   each coupled state, and orthogonality b/w different J's at same M.
*/

#include "../core/complex.h"
#include "../core/vector.h"
#include "../physics/angular.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.6f expected=%.6f err=%.2e\n", label, got, expected, err);
  return err > tol;
}

// Two spin-1/2 -> triplet/singlet
static int test_two_spin_half(void) {
  int fail = 0;
  double tol = 1e-12;
  double inv_sqrt2 = 1.0 / sqrt(2.0);

  // Basis index convention: i1*2+i2, i=0 -> m=-1/2 (down), i=1 -> m=+1/2 (up)
  // Triplet M=+1: |up up> = 1
  cvector_t *v = couple_states(1, 1, 2, 2);
  if (!v) {
    printf("  FAIL: couple_states(1,1,2,2) returned NULL\n");
    return 1;
  }
  fail |= check_close(v->data[3].re, 1.0, tol, "triplet M=+1, |up up>");
  fail |= check_close(v->data[0].re, 0.0, tol, "triplet M=+1, |down down>");
  cvector_free(v);

  // Triplet M=0: (|up down> + |down up>)/sqrt2
  v = couple_states(1, 1, 2, 0);
  if (!v) {
    printf("  FAIL: couple_states(1,1,2,0) returned NULL\n");
    return 1;
  }
  fail |= check_close(v->data[1].re, inv_sqrt2, tol, "triplet M=0, |down up>");
  fail |= check_close(v->data[2].re, inv_sqrt2, tol, "triplet M=0, |up down>");
  cvector_free(v);

  // Singlet M=0: (|up down> - |down up>)/sqrt2
  v = couple_states(1, 1, 0, 0);
  if (!v) {
    printf("  FAIL: couple_states(1,1,0,0) returned NULL\n");
    return 1;
  }
  fail |= check_close(v->data[1].re, -inv_sqrt2, tol, "singlet, |down up>");
  fail |= check_close(v->data[2].re, inv_sqrt2, tol, "singlet, |up down>");
  cvector_free(v);

  // Triplet M=-1: |down down> = 1
  v = couple_states(1, 1, 2, -2);
  if (!v) {
    printf("  FAIL: couple_states(1,1,2,-2) returned NULL\n");
    return 1;
  }
  fail |= check_close(v->data[0].re, 1.0, tol, "triplet M=-1, |down down>");
  cvector_free(v);

  return fail;
}

// j1=1 (orbital l=1), j2=1/2 (spin) -> J=3/2 or J=1/2.
static int test_l1_s_half(void) {
  int fail = 0;
  double tol = 1e-10;

  int J2_list[8];
  int count = couple_allowed_J(2, 1, J2_list);
  printf("  allowed J (doubled) for j1=1,j2=1/2: count=%d ->", count);
  for (int i = 0; i < count; i++)
    printf(" %d", J2_list[i]);
  printf("\n");
  if (count != 2 || J2_list[0] != 1 || J2_list[1] != 3) {
    printf("  FAIL: expected J_2 in {1,3}\n");
    fail = 1;
  }

  // Normalization: sum of squares of a coupled state must be 1.
  int M2 = 1; // M=1/2, valid for both J=1/2 and J=3/2
  for (int idx = 0; idx < count; idx++) {
    int J2 = J2_list[idx];
    cvector_t *v = couple_states(2, 1, J2, M2);
    if (!v) {
      printf("  FAIL: couple_states(2,1,%d,%d) returned NULL\n", J2, M2);
      fail = 1;
      continue;
    }
    double norm_sq = 0.0;
    for (int i = 0; i < v->n; i++)
      norm_sq += v->data[i].re * v->data[i].re;
    fail |= check_close(norm_sq, 1.0, tol, "normalization");
    cvector_free(v);
  }

  // Orthogonality: two J states at same M must be orthogonal.
  cvector_t *v1 = couple_states(2, 1, J2_list[0], M2);
  cvector_t *v2 = couple_states(2, 1, J2_list[1], M2);
  if (v1 && v2 && v1->n == v2->n) {
    double dot = 0.0;
    for (int i = 0; i < v1->n; i++)
      dot += v1->data[i].re * v2->data[i].re;
    fail |= check_close(dot, 0.0, tol, "orthogonality J=1/2 vs J=3/2");
  } else {
    fail = 1;
  }
  cvector_free(v1);
  cvector_free(v2);

  return fail;
}

int main(void) {
  int failed = 0;

  printf("Two spin-1/2 -> triplet/singlet:\n");
  failed += test_two_spin_half();

  printf("l=1 (x) s=1/2 (spin-orbit basis):\n");
  failed += test_l1_s_half();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");
  return 0;
}
