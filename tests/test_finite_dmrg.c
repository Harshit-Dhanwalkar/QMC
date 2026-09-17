/*
 * Test: finite-system DMRG sweeps for the open-boundary spin-1/2 XXZ chain
 *
 * NOTE: Validates against an independent brute-force exact diagonalization,
 * checks that repeated sweeps never make reported central-bond energy worse
 * (variational principle: every sweep step is itself an exact ground-state
 * diagonalization within a subspace, so refining that subspace can only lower
 * or hold energy), and spot-checks that finite sweeping matches or improves on
 * infinite-system algorithm at small, truncation-limited m
 */

#include "../core/complex.h"
#include "../physics/dmrg.h"
#include "../physics/finite_dmrg.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef RUNNING_ON_VALGRIND
#define RUNNING_ON_VALGRIND 0
#endif

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);
  if (err > tol) {
    printf("  FAIL: %s (got=%.10f expected=%.10f err=%.2e > tol=%.2e)\n", label,
           got, expected, err, tol);
    failures++;
  } else {
    printf("  ok:   %s (err=%.2e)\n", label, err);
  }
}

static void check(int cond, const char *label) {
  if (!cond) {
    printf("  FAIL: %s\n", label);
    failures++;
  } else {
    printf("  ok:   %s\n", label);
  }
}

/* Brute-force OBC XXZ ground energy, independent of both dmrg.c and
 * finite_dmrg.c. Only tractable for small N; used here for N<=12. */
static double brute_force_ed(int N, double Jz, double Jxy) {
  cmatrix_t *Sz = cmatrix_alloc(2, 2);
  CMAT(Sz, 0, 0) = c_real(0.5);
  CMAT(Sz, 0, 1) = c_zero();
  CMAT(Sz, 1, 0) = c_zero();
  CMAT(Sz, 1, 1) = c_real(-0.5);

  cmatrix_t *Sp = cmatrix_alloc(2, 2);
  CMAT(Sp, 0, 0) = c_zero();
  CMAT(Sp, 0, 1) = c_real(1.0);
  CMAT(Sp, 1, 0) = c_zero();
  CMAT(Sp, 1, 1) = c_zero();

  cmatrix_t *Sm = cmatrix_adjoint(Sp);
  cmatrix_t *I2 = cmatrix_alloc(2, 2);

  CMAT(I2, 0, 0) = c_real(1.0);
  CMAT(I2, 0, 1) = c_zero();
  CMAT(I2, 1, 0) = c_zero();
  CMAT(I2, 1, 1) = c_real(1.0);

  int dim = 1 << N;
  cmatrix_t *H = cmatrix_alloc(dim, dim);
  for (int i = 0; i < dim * dim; i++) {
    H->data[i] = c_zero();
  }

  for (int j = 0; j < N - 1; j++) {
    for (int which = 0; which < 3; which++) {
      cmatrix_t *acc = NULL;

      for (int s = 0; s < N; s++) {
        const cmatrix_t *op;

        if (s == j) {
          op = (which == 0) ? Sz : (which == 1) ? Sp : Sm;
        } else if (s == j + 1) {
          op = (which == 0) ? Sz : (which == 1) ? Sm : Sp;
        } else {
          op = I2;
        }

        if (!acc) {
          acc = cmatrix_copy(op);
        } else {
          cmatrix_t *nacc = cmatrix_kron(acc, op);

          cmatrix_free(acc);
          acc = nacc;
        }
      }

      double coeff = (which == 0) ? Jz : (Jxy / 2.0);
      for (int i = 0; i < dim * dim; i++) {
        H->data[i] = c_add(H->data[i], c_scale(acc->data[i], coeff));
      }

      cmatrix_free(acc);
    }
  }

  eigen_t *eig = cmatrix_eigh(H);
  double E0 = eig->eigenvalues[0];

  eigen_free(eig);
  cmatrix_free(H);
  cmatrix_free(Sz);
  cmatrix_free(Sp);
  cmatrix_free(Sm);
  cmatrix_free(I2);

  return E0;
}

/*
 * Cross-check against brute-force ED for small N (N=8), isotropic
 * Heisenberg AFM (Jz=Jxy=1): energy error should shrink as m grows, and
 * reach near machine precision once m is no longer restrictive
 */
static void test_vs_brute_force_ed(void) {
  printf("  === Finite DMRG vs brute-force exact diagonalization (N=8) ===\n");
  double Jz = 1.0, Jxy = 1.0;
  int N = 8;
  double E_exact = brute_force_ed(N, Jz, Jxy);

  /* NOTE: m=8's tolerance is looser than test_dmrg.c's equivalent case: for
   * N=8 the exact central-cut bond dimension is 2^(N/2)=16, so m=8 still
   * discards real information here and cannot be lossless the way
   * dmrg_run(8, ..., 8) incidentally is (that function returns its
   * superblock energy one step before ever applying an m=8 truncation, so
   * it never actually pays the cost its own requested m implies for this
   * particular N; finite_dmrg_run intentionally truncates every block,
   * including the last, so it does pay that cost here). */
  struct {
    int m;
    double tol;
  } ms[] = {{2, 2.5}, {4, 0.1}, {8, 5e-6}, {16, 1e-6}};

  double prev_err = 1e300;
  for (size_t i = 0; i < sizeof(ms) / sizeof(ms[0]); i++) {
    int m = ms[i].m;
    finite_dmrg_result_t *r = finite_dmrg_run(N, Jz, Jxy, m, 4);
    check(r != NULL, "finite_dmrg_run succeeds");
    if (!r) {
      continue;
    }

    check(r->N_reached == N, "N_reached == N_target (N even)");
    double err = fabs(r->energy - E_exact);
    char label[128];

    snprintf(label, sizeof label, "N=%d m=%d matches brute-force ED", N, m);

    check_close(r->energy, E_exact, ms[i].tol, label);
    check(err <= prev_err + 1e-9, "energy error is non-increasing as m grows");

    prev_err = err;

    finite_dmrg_result_free(r);
  }
}

/*
 * Reference values independently cross-checked (same values test_dmrg.c
 * uses for infinite algorithm): finite sweeps at generous m must reach same
 * lossless answer
 */
static void test_reference_values(void) {
  printf("  === Finite DMRG vs validated reference energies (N=10,12) ===\n");
  double Jz = 1.0;
  double Jxy = 1.0;

  struct {
    int N;
    double E_exact;
  } cases[] = {
      {10, -4.25803521},
      {12, -5.14209063},
  };

  for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
    int N = cases[c].N;
    int m = 24;
    finite_dmrg_result_t *r = finite_dmrg_run(N, Jz, Jxy, m, 4);
    check(r != NULL, "finite_dmrg_run succeeds");
    if (!r) {
      continue;
    }

    char label[128];
    snprintf(label, sizeof label, "N=%d lossless-m matches reference", N);

    check_close(r->energy, cases[c].E_exact, 1e-6, label);
    check(r->truncation_error < 1e-8, "truncation_error ~0 at lossless m");

    finite_dmrg_result_free(r);
  }
}

/*
 * Every sweep step is an exact ground-state diagonalization of current
 * superblock, so the central-bond checkpoint energy recorded after each
 * successive sweep must never increase (up to solver round-off)
 */
static void test_sweep_energy_monotonic(void) {
  printf("  === Central-bond energy is non-increasing across sweeps (N=14) "
         "===\n");
  double Jz = 1.0, Jxy = 1.0;
  int N = 14;
  int m = RUNNING_ON_VALGRIND ? 4 : 6;
  int n_sweeps = RUNNING_ON_VALGRIND ? 2 : 5;

  finite_dmrg_result_t *r = finite_dmrg_run(N, Jz, Jxy, m, n_sweeps);
  check(r != NULL, "finite_dmrg_run succeeds");
  if (!r) {
    return;
  }

  check(r->n_sweeps == n_sweeps, "n_sweeps reported matches request");
  check(r->sweep_energy != NULL, "sweep_energy array populated");

  if (r->sweep_energy) {
    check_close(r->sweep_energy[n_sweeps - 1], r->energy, 1e-14,
                "sweep_energy[last] == reported energy");

    for (int s = 1; s < n_sweeps; s++) {
      char label[128];
      snprintf(label, sizeof label, "sweep %d energy <= sweep %d energy", s,
               s - 1);
      check(r->sweep_energy[s] <= r->sweep_energy[s - 1] + 1e-9, label);
    }
  }

  finite_dmrg_result_free(r);
}

/*
 * At small, truncation-limited m, finite sweeping should match or improve on
 * infinite algorithm's result for same N and m: sweeping only ever re-optimizes
 * previously-built blocks using better information, whereas infinite algorithm
 * builds every block once and never revisits it. Uses a modest tolerance since
 * the two algorithms' final steps are not directly nested subspaces
 */
static void test_matches_or_improves_infinite(void) {
  printf("  === Finite sweeps match/improve on infinite algorithm (small m) "
         "===\n");
  double Jz = 1.0, Jxy = 1.0;
  int N = 16;
  int m = 6;

  dmrg_result_t *r_inf = dmrg_run(N, Jz, Jxy, m);
  finite_dmrg_result_t *r_fin = finite_dmrg_run(N, Jz, Jxy, m, 4);

  check(r_inf != NULL, "dmrg_run succeeds");
  check(r_fin != NULL, "finite_dmrg_run succeeds");

  if (r_inf && r_fin) {
    check(r_fin->energy <= r_inf->energy + 1e-6,
          "finite energy <= infinite energy + tol");
  }

  free(r_inf);
  finite_dmrg_result_free(r_fin);
}

static void test_invalid_input(void) {
  printf(" === Invalid input handling ===\n");
  check(finite_dmrg_run(2, 1.0, 1.0, 8, 2) == NULL, "N_target < 4 rejected");
  check(finite_dmrg_run(8, 1.0, 1.0, 0, 2) == NULL, "m_max < 1 rejected");
  check(finite_dmrg_run(8, 1.0, 1.0, 8, 0) == NULL, "n_sweeps < 1 rejected");
  check(finite_dmrg_run(8, NAN, 1.0, 8, 2) == NULL, "non-finite Jz rejected");
  check(finite_dmrg_run(8, 1.0, INFINITY, 8, 2) == NULL,
        "non-finite Jxy rejected");
}

int main(void) {
  printf(" > Finite-system DMRG: XXZ open chain sweep tests\n");

  test_vs_brute_force_ed();
  test_reference_values();
  test_sweep_energy_monotonic();
  test_matches_or_improves_infinite();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_finite_dmrg checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_finite_dmrg check(s) FAILED.\n", failures);
    return 1;
  }
}
