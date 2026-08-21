/*
 * Property-based / fuzz testing for physics/molecular_integrals.c.
 *
 * NOTE: Rather than comparing against fixed reference numbers, this generates
 * many random basis functions/geometries and checks mathematical invariants
 * that must hold for any valid input -> like symmetry, positivity, translation
 * invariance, and known closed-form  special cases.
 *
 * Uses a fixed PRNG seed for reproducibility (a failing fuzz run should be
 * reproducible) but draws many (n_trials) random cases per property, so this
 * still explores far more of the input space than any fixed hand-written test
 * geometry.
 */

#include "../core/matrix.h"
#include "../core/random.h"
#include "../physics/molecular_integrals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;
static int checks = 0;

static void check_true(int cond, const char *label) {
  checks++;
  if (!cond) {
    printf("  FAIL: %s\n", label);
    failures++;
  }
}

static void check_close(double got, double expected, double tol,
                        const char *label) {
  checks++;
  double err = fabs(got - expected);
  if (err > tol) {
    printf("  FAIL: %s: got=%.12e expected=%.12e err=%.3e (tol=%.3e)\n", label,
           got, expected, err, tol);
    failures++;
  }
}

/* ----- random basis function generation ----- */

static void random_center(rng_state_t *rng, double c[3], double spread) {
  for (int d = 0; d < 3; d++) {
    c[d] = rng_uniform_range(rng, -spread, spread);
  }
}

/*
 * NOTE: Random contracted s/p/d basis function: 1-3 primitives, exponents in a
 * physically-reasonable range (0.1 to 10, covering both diffuse and fairly
 * tight primitives : not pushing into x>250 core-orbital regime here, since
 * that needs same-center pairs specifically, handled separately below), random
 * (but not all-zero) contraction coefficients.
 */
static basis_function_t *
random_basis_function(rng_state_t *rng, const double center[3], int max_l) {
  int l = (int)rng_uniform_range(rng, 0, max_l + 1);
  int m = (int)rng_uniform_range(rng, 0, max_l - l + 1);
  int n = (int)rng_uniform_range(rng, 0, max_l - l - m + 1);
  if (l + m + n > max_l) {
    n = max_l - l - m;

    if (n < 0) {
      n = 0;
    }
  }

  int n_prim = 1 + (int)rng_uniform_range(rng, 0, 3);
  double exps[3], coefs[3];
  for (int p = 0; p < n_prim; p++) {
    exps[p] = pow(10.0, rng_uniform_range(rng, -1.0, 1.0)); // 0.1 .. 10
    coefs[p] = rng_uniform_range(rng, 0.2, 1.0);
  }

  basis_function_t *bf =
      basis_function_alloc(l, m, n, center, n_prim, exps, coefs);
  molint_normalize_contraction(bf);

  return bf;
}

/* ----- Boys function properties ----- */

static void test_boys_function_properties(void) {
  printf("test_boys_function_properties:\n");
  rng_state_t rng;
  rng_seed(&rng, 12345);

  // F_n(0) = 1/(2n+1) exactly
  for (int n = 0; n <= 8; n++) {
    check_close(boys_function(n, 0.0), 1.0 / (2 * n + 1), 1e-12,
                "F_n(0) = 1/(2n+1)");
  }

  /* Monotonically decreasing in x for fixed n (F_n(x) is an integral of a
   * positive integrand times a monotonically-decreasing-in-x weight) */
  for (int trial = 0; trial < 200; trial++) {
    int n = (int)rng_uniform_range(&rng, 0, 6);
    double x1 = rng_uniform_range(&rng, 0.0, 500.0);
    double x2 = x1 + rng_uniform_range(&rng, 0.01, 50.0);
    double f1 = boys_function(n, x1);
    double f2 = boys_function(n, x2);

    check_true(f2 <= f1 + 1e-13, "F_n(x) monotonically non-increasing in x");
  }

  /* Positivity for every n, x */
  for (int trial = 0; trial < 200; trial++) {
    int n = (int)rng_uniform_range(&rng, 0, 8);
    double x = rng_uniform_range(&rng, 0.0, 600.0);

    check_true(boys_function(n, x) >= 0.0, "F_n(x) >= 0");
  }

  /*
   * NOTE: Continuity across the x=30 series/closed-form switching boundary.
   * Compare boys_function_array's F_0..F_8 just below and just above the switch
   * : should agree to near machine precision if two regimes are correctly
   * stitched together.
   */
  double F_below[9], F_above[9];

  boys_function_array(8, 29.999, F_below);
  boys_function_array(8, 30.001, F_above);

  for (int n = 0; n <= 8; n++) {
    char label[64];
    snprintf(label, sizeof label, "F_%d continuous across x=30 boundary", n);

    check_close(F_above[n], F_below[n], 1e-4, label);
  }

  /*
   * NOTE: Recursion consistency: boys_function(n,x) (presumably an independent
   * single-value path or a thin wrapper) should agree with
   * boys_function_array's batch computation at the same n, x, across both
   * regimes (n small, x small and n moderate, x large past the switch which is
   * irectly exercising the large-x).
   */
  for (int trial = 0; trial < 100; trial++) {
    int nmax = (int)rng_uniform_range(&rng, 1, 8);
    double x = rng_uniform_range(&rng, 0.0, 600.0); /* spans both regimes */
    double *F = malloc((size_t)(nmax + 1) * sizeof(double));

    boys_function_array(nmax, x, F);

    for (int n = 0; n <= nmax; n++) {
      check_close(boys_function(n, x), F[n], 1e-8,
                  "boys_function matches boys_function_array at same (n,x)");
    }

    free(F);
  }
}

/* ----- symmetry / positivity / translation invariance ----- */

static void test_overlap_properties(void) {
  printf("test_overlap_properties:\n");
  rng_state_t rng;
  rng_seed(&rng, 2024);

  for (int trial = 0; trial < 100; trial++) {
    double ca[3], cb[3];

    random_center(&rng, ca, 3.0);
    random_center(&rng, cb, 3.0);

    basis_function_t *a = random_basis_function(&rng, ca, 2);
    basis_function_t *b = random_basis_function(&rng, cb, 2);

    double Sab = gto_overlap(a, b);
    double Sba = gto_overlap(b, a);
    check_close(Sab, Sba, 1e-10, "S(a,b) = S(b,a)");

    double Saa = gto_overlap(a, a);
    check_true(Saa > 0.0, "S(a,a) > 0");

    basis_function_free(a);
    basis_function_free(b);
  }

  /* normalized self-overlap should be exactly 1 */
  for (int trial = 0; trial < 50; trial++) {
    double c[3];
    random_center(&rng, c, 3.0);

    basis_function_t *a = random_basis_function(&rng, c, 2);
    check_close(gto_overlap(a, a), 1.0, 1e-8,
                "normalized self-overlap S(a,a) = 1");

    basis_function_free(a);
  }
}

static void test_kinetic_properties(void) {
  printf("test_kinetic_properties:\n");
  rng_state_t rng;
  rng_seed(&rng, 999);

  for (int trial = 0; trial < 100; trial++) {
    double ca[3], cb[3];

    random_center(&rng, ca, 3.0);
    random_center(&rng, cb, 3.0);

    basis_function_t *a = random_basis_function(&rng, ca, 2);
    basis_function_t *b = random_basis_function(&rng, cb, 2);

    double Tab = gto_kinetic(a, b);
    double Tba = gto_kinetic(b, a);
    check_close(Tab, Tba, 1e-9, "T(a,b) = T(b,a)");

    double Taa = gto_kinetic(a, a);
    check_true(Taa >= 0.0, "T(a,a) >= 0 (kinetic energy expectation)");

    basis_function_free(a);
    basis_function_free(b);
  }
}

static void test_eri_symmetry(void) {
  printf("test_eri_symmetry:\n");
  rng_state_t rng;
  rng_seed(&rng, 42);

  for (int trial = 0; trial < 60; trial++) {
    double c1[3], c2[3], c3[3], c4[3];

    random_center(&rng, c1, 2.0);
    random_center(&rng, c2, 2.0);
    random_center(&rng, c3, 2.0);
    random_center(&rng, c4, 2.0);

    basis_function_t *a = random_basis_function(&rng, c1, 1);
    basis_function_t *b = random_basis_function(&rng, c2, 1);
    basis_function_t *c = random_basis_function(&rng, c3, 1);
    basis_function_t *d = random_basis_function(&rng, c4, 1);

    double abcd = gto_eri(a, b, c, d);
    double bacd = gto_eri(b, a, c, d);
    double abdc = gto_eri(a, b, d, c);
    double cdab = gto_eri(c, d, a, b);
    double dcba = gto_eri(d, c, b, a);

    double tol = 1e-8 * (1.0 + fabs(abcd));

    check_close(bacd, abcd, tol, "(ab|cd) = (ba|cd)");
    check_close(abdc, abcd, tol, "(ab|cd) = (ab|dc)");
    check_close(cdab, abcd, tol, "(ab|cd) = (cd|ab)");
    check_close(dcba, abcd, tol, "(ab|cd) = (dc|ba)");

    basis_function_free(a);
    basis_function_free(b);
    basis_function_free(c);
    basis_function_free(d);
  }

  // Coulomb self-repulsion (aa|aa) must be strictly positive
  for (int trial = 0; trial < 30; trial++) {
    double c[3];
    random_center(&rng, c, 2.0);

    basis_function_t *a = random_basis_function(&rng, c, 2);
    double aaaa = gto_eri(a, a, a, a);

    check_true(aaaa > 0.0, "(aa|aa) > 0 (self-repulsion)");
    basis_function_free(a);
  }
}

static void test_translation_invariance(void) {
  printf("test_translation_invariance:\n");
  rng_state_t rng;
  rng_seed(&rng, 7);

  for (int trial = 0; trial < 40; trial++) {
    double ca[3], cb[3], cC[3], shift[3];
    random_center(&rng, ca, 2.0);
    random_center(&rng, cb, 2.0);
    random_center(&rng, cC, 2.0);
    random_center(&rng, shift, 5.0);

    basis_function_t *a = random_basis_function(&rng, ca, 1);
    basis_function_t *b = random_basis_function(&rng, cb, 1);

    double ca_shift[3], cb_shift[3], cC_shift[3];
    for (int d = 0; d < 3; d++) {
      ca_shift[d] = ca[d] + shift[d];
      cb_shift[d] = cb[d] + shift[d];
      cC_shift[d] = cC[d] + shift[d];
    }
    basis_function_t *a_shift =
        basis_function_alloc(a->l, a->m, a->n, ca_shift, a->n_primitives,
                             a->exponents, a->coefficients);
    basis_function_t *b_shift =
        basis_function_alloc(b->l, b->m, b->n, cb_shift, b->n_primitives,
                             b->exponents, b->coefficients);

    check_close(gto_overlap(a, b), gto_overlap(a_shift, b_shift), 1e-9,
                "overlap is translation invariant");
    check_close(gto_kinetic(a, b), gto_kinetic(a_shift, b_shift), 1e-8,
                "kinetic is translation invariant");
    check_close(gto_nuclear_attraction(a, b, cC),
                gto_nuclear_attraction(a_shift, b_shift, cC_shift), 1e-8,
                "nuclear attraction is translation invariant "
                "(nuclear center shifted along with the basis functions)");

    basis_function_free(a);
    basis_function_free(b);
    basis_function_free(a_shift);
    basis_function_free(b_shift);
  }
}

/* ----- Independent numerical-quadrature cross-check for overlap ----- */

/* NOTE: Brute-force 3D trapezoidal-rule numerical overlap integral, completely
 * independent of the analytic McMurchie-Davidson implementation being tested :
 * deliberately simple/slow/low-precision (coarse grid), just enough to catch a
 * gross sign or factor-of-N error the symmetry checks above couldn't (those are
 * self-consistency checks on gto_overlap alone). */
static double numerical_overlap(const basis_function_t *a,
                                const basis_function_t *b, double half_extent,
                                int n_grid) {
  double total = 0.0;
  double h = 2.0 * half_extent / n_grid;

  for (int ix = 0; ix < n_grid; ix++) {
    double x = -half_extent + (ix + 0.5) * h;

    for (int iy = 0; iy < n_grid; iy++) {
      double y = -half_extent + (iy + 0.5) * h;

      for (int iz = 0; iz < n_grid; iz++) {
        double z = -half_extent + (iz + 0.5) * h;

        double va = 0.0, vb = 0.0;
        double dxa = x - a->center[0], dya = y - a->center[1],
               dza = z - a->center[2];
        double ra2 = dxa * dxa + dya * dya + dza * dza;
        double poly_a = pow(dxa, a->l) * pow(dya, a->m) * pow(dza, a->n);

        for (int p = 0; p < a->n_primitives; p++) {
          va += a->coefficients[p] * exp(-a->exponents[p] * ra2);
        }

        va *= poly_a;

        double dxb = x - b->center[0], dyb = y - b->center[1],
               dzb = z - b->center[2];
        double rb2 = dxb * dxb + dyb * dyb + dzb * dzb;
        double poly_b = pow(dxb, b->l) * pow(dyb, b->m) * pow(dzb, b->n);

        for (int p = 0; p < b->n_primitives; p++) {
          vb += b->coefficients[p] * exp(-b->exponents[p] * rb2);
        }

        vb *= poly_b;

        total += va * vb;
      }
    }
  }

  return total * h * h * h;
}

static void test_overlap_vs_numerical_quadrature(void) {
  printf("test_overlap_vs_numerical_quadrature:\n");
  rng_state_t rng;
  rng_seed(&rng, 555);

  /* Small number of trials, kept coarse grid and only s/p functions close
   * together so grid resolves Gaussians without needing an large n_grid
   */
  for (int trial = 0; trial < 6; trial++) {
    double ca[3] = {0, 0, 0};
    double cb[3];
    random_center(&rng, cb, 1.0);

    basis_function_t *a = random_basis_function(&rng, ca, 1);
    basis_function_t *b = random_basis_function(&rng, cb, 1);

    double analytic = gto_overlap(a, b);
    double numerical = numerical_overlap(a, b, 6.0, 80);

    /* coarse grid -> loose tolerance, this is a sanity cross-check not a
     * precision benchmark
     */
    check_close(analytic, numerical, 5e-3,
                "analytic overlap matches independent numerical quadrature");

    basis_function_free(a);
    basis_function_free(b);
  }
}

int main(void) {
  test_boys_function_properties();
  test_overlap_properties();
  test_kinetic_properties();
  test_eri_symmetry();
  test_translation_invariance();
  test_overlap_vs_numerical_quadrature();

  printf("\n%d checks run, %d failed\n", checks, failures);
  if (failures == 0) {
    printf("All molecular_integrals fuzz checks PASSED\n");
    return 0;
  }
  return 1;
}
