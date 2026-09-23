/*
 * Test: physics/spin_chain.c's spin_dsf_continued_fraction (Lanczos
 * continued-fraction evaluation of dynamical structure factor
 * S^{zz}(q, \omega))
 *
 * test_spin_chain.c already has one DSF check (N=6, momentum-transfer index
 * q_index=1, sum rule vs I0). This file covers three properties that check
 * doesn't:
 *
 *   1. The q_index=0 special case: for a state in S^z_total=0 sector (which
 *      Heisenberg antiferromagnet's true ground state always is, Lieb-Mattis),
 *      S^z_{q=0} is proportional to S^z_total operator, which annihilates it
 *      exactly. So I0 (and hence whole DSF) must come out at machine-precision
 *      zero here - an exact, not just approximate, sum-rule value, on a code
 *      path (q_index=0 => target sector == source sector) q_index=1 test never
 *      exercises
 *   2. Positivity: S(q, \omega) is (up to broadening) a sum of positive
 *      spectral weights times a positive Lorentzian, so it must be non-negative
 *      at every frequency. A sign error anywhere in continued-fraction
 *      recursion could still pass a sum-rule check (total area could come out
 *      right by cancellation) while violating this pointwise
 *   3. An independent (N, q) sum-rule check on a second, smaller ring (N=4)
 *      than existing N=6 test, as a cross-check that sum rule isn't an accident
 *      of one particular system size
 */

#include "../core/sparse.h"
#include "../core/vector.h"
#include "../physics/spin_chain.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check(int cond, const char *msg) {
  printf("  %s: %s\n", msg, cond ? "ok" : "FAILED");
  if (!cond) {
    failures++;
  }
}

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  double err = fabs(got - expected);
  printf("  %s: got=%.10f expected=%.10f err=%.2e\n", msg, got, expected, err);
  if (err > tol) {
    failures++;
  }
}

/* Scan every k in the nup=N/2 (S^z_total=0) sector and return the momentum
 * index of the lowest-energy state - the Heisenberg ring's true ground
 * state always sits in this sector (Lieb-Mattis), so this is enough to find
 * it without also scanning nup */
static int find_ground_k_at_half_filling(int N, double *E0_out) {
  int nup = N / 2;
  double best = 1e300;
  int best_k = -1;

  for (int k = 0; k < N; k++) {
    spin_sector_t *sec = spin_sector_build(N, nup, k);
    if (!sec || sec->dim == 0) {
      spin_sector_free(sec);
      continue;
    }

    sparse_matrix_t *H = spin_sector_hamiltonian(sec, 1.0, 1.0, 1);
    lanczos_result_t *res =
        lanczos_eigs(H, 1, sec->dim < 60 ? sec->dim : 60, 1e-12);

    if (res && res->values[0] < best) {
      best = res->values[0];
      best_k = k;
    }

    lanczos_free(res);
    sparse_free(H);
    spin_sector_free(sec);
  }

  if (E0_out) {
    *E0_out = best;
  }
  return best_k;
}

/* Runs the full pipeline (ground state -> S^z_q|psi0> -> Lanczos
 * tridiagonalization -> continued-fraction DSF) for a given N/q_index, and
 * hands back E0, I0, and the tridiagonal coefficients for the caller to
 * evaluate spin_dsf_continued_fraction with. Returns 0 on success */
static int build_dsf_pipeline(int N, int q_index, double *E0_out,
                              double *I0_out, lanczos_tridiag_t **tri_out) {
  double E0;
  int k_gs = find_ground_k_at_half_filling(N, &E0);
  if (k_gs < 0) {
    return 1;
  }

  int nup = N / 2;
  spin_sector_t *sec0 = spin_sector_build(N, nup, k_gs);
  sparse_matrix_t *H0 = spin_sector_hamiltonian(sec0, 1.0, 1.0, 1);
  lanczos_result_t *gs = lanczos_eigs(H0, 1, sec0->dim, 1e-12);

  cvector_t *psi0 = cvector_alloc(sec0->dim);
  for (int i = 0; i < sec0->dim; i++) {
    psi0->data[i] = CMAT(gs->vectors, i, 0);
  }

  spin_sector_t *target = NULL;
  cvector_t *phi0 = NULL;
  double I0 = -1.0;
  int rc = spin_apply_szq(sec0, psi0, q_index, &target, &phi0, &I0);

  int status = 1;
  lanczos_tridiag_t *tri = NULL;

  if (rc == 0) {
    if (I0 > 1e-14) {
      cvector_t *f0 = cvector_copy(phi0);
      cvector_normalize(f0);

      sparse_matrix_t *Ht = spin_sector_hamiltonian(target, 1.0, 1.0, 1);
      tri = lanczos_tridiagonalize(Ht, f0, target->dim, 1e-12);
      sparse_free(Ht);
      cvector_free(f0);

      status = (tri != NULL) ? 0 : 1;
    } else {
      /* I0 essentially zero (e.g. q_index=0): nothing to tridiagonalize,
       * the caller only needs I0 itself */
      status = 0;
    }
  }

  *E0_out = gs->values[0];
  *I0_out = I0;
  *tri_out = tri;

  if (phi0) {
    cvector_free(phi0);
  }
  if (target) {
    spin_sector_free(target);
  }
  cvector_free(psi0);
  lanczos_free(gs);
  sparse_free(H0);
  spin_sector_free(sec0);

  return status;
}

static void test_q_zero_annihilates_ground_state(void) {
  printf("  === Test q zero annihilates ground state ===\n");
  printf("  (S^z_q=0 is proportional to S^z_total, which is exactly zero on "
         "the S^z_total=0 ground state - I0 must vanish at machine "
         "precision)\n");

  int N = 6;
  double E0;
  double I0;
  lanczos_tridiag_t *tri = NULL;
  int rc = build_dsf_pipeline(N, 0, &E0, &I0, &tri);

  check(rc == 0, "pipeline runs without error at q_index=0");
  check_close(I0, 0.0, 1e-10, "I0 is (machine-precision) zero at q_index=0");

  // With I0 essentially zero, the DSF itself must be zero everywhere too -
  // there's no spectral weight left for any (alpha, beta) recursion to
  // redistribute. tri is NULL here (build_dsf_pipeline skips
  // tridiagonalization when I0 ~ 0), so evaluate directly with m=0
  double S = spin_dsf_continued_fraction(NULL, NULL, 0, E0, I0, 1.0, 0.1);
  check_close(S, 0.0, 1e-12, "S(q=0, \\omega) is exactly zero (I0=0, m=0)");

  if (tri) {
    lanczos_tridiag_free(tri);
  }
}

static void test_dsf_is_nonnegative(void) {
  printf("  === Test DSF is nonnegative ===\n");
  printf("  (S(q, \\omega) is a sum of positive Lorentzians weighted by "
         "positive spectral weights - must never go negative)\n");

  int N = 6;
  int q_index = 1;
  double E0, I0;
  lanczos_tridiag_t *tri = NULL;
  int rc = build_dsf_pipeline(N, q_index, &E0, &I0, &tri);

  check(rc == 0, "pipeline runs without error");
  check(tri != NULL, "Lanczos tridiagonalization succeeded");

  if (tri) {
    double eta = 0.05;
    double omega_min = -4.0;
    double omega_max = 10.0;
    int npts = 4000;
    double min_S = 1e300;

    for (int i = 0; i <= npts; i++) {
      double w = omega_min + i * (omega_max - omega_min) / npts;
      double S = spin_dsf_continued_fraction(tri->alpha, tri->beta, tri->m, E0,
                                             I0, w, eta);
      if (S < min_S) {
        min_S = S;
      }
    }

    printf("  min S(q, \\omega) over [%.1f, %.1f]: %.3e\n", omega_min,
           omega_max, min_S);
    // Small negative tolerance for floating-point noise only
    check(min_S > -1e-9, "S(q, \\omega) is non-negative across the sampled "
                         "frequency range");

    lanczos_tridiag_free(tri);
  }
}

static void test_n4_sum_rule_cross_check(void) {
  printf("  === Test N4 sum rule cross check ===\n");
  printf("  (independent sum-rule check on a second, smaller ring than "
         "test_spin_chain.c's N=6 case)\n");

  int N = 4;
  int q_index = 1;
  double E0, I0;
  lanczos_tridiag_t *tri = NULL;
  int rc = build_dsf_pipeline(N, q_index, &E0, &I0, &tri);

  check(rc == 0, "pipeline runs without error");
  check(I0 > 0.0,
        "I0 is strictly positive for q_index=1 (nontrivial spectral weight)");
  check(tri != NULL, "Lanczos tridiagonalization succeeded");

  if (tri) {
    double eta = 0.05;
    double omega_min = -2.0;
    double omega_max = 8.0;
    int npts = 20000;
    double domega = (omega_max - omega_min) / npts;
    double integral = 0.0;

    for (int i = 0; i < npts; i++) {
      double w = omega_min + (i + 0.5) * domega;
      double S = spin_dsf_continued_fraction(tri->alpha, tri->beta, tri->m, E0,
                                             I0, w, eta);
      integral += S * domega;
    }

    check_close(
        integral, I0, 5e-3,
        "numerical integral of S(q,omega) over \\omega matches I0  (N=4)");

    lanczos_tridiag_free(tri);
  }
}

int main(void) {
  test_q_zero_annihilates_ground_state();
  test_dsf_is_nonnegative();
  test_n4_sum_rule_cross_check();

  if (failures > 0) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll checks passed\n");

  return 0;
}
