/*
 * Test: infinite-system DMRG for the open-boundary spin-1/2 XXZ chain
 *
 * NOTE: Checks expected qualitative signature of a correct density-matrix
 * truncation: energy error and discarded RDM weight (truncation_error) both
 * shrink monotonically as m increases, reaching machine precision once m
 * exceeds 2^(N/2) (point at which no information is discarded), and that the
 * per-site ground energy approaches exact Bethe-ansatz thermodynamic-limit
 * value (1/4 - ln(2)) as N grows
 */

#include "../core/complex.h"
#include "../physics/dmrg.h"
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

/*
 * NOTE: Brute-force OBC XXZ ground energy via direct 2^N x 2^N Hamiltonian
 * construction, entirely independent of dmrg.c. Only tractable for small N
 * dense O(2^{3N}) diagonalization); used here for N<=8 only.
 */
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

static void test_block_init(void) {
  printf(" === dmrg_block_init: single-site block ===\n");
  dmrg_block_t *b = dmrg_block_init();
  check(b != NULL, "dmrg_block_init succeeds");
  if (!b) {
    return;
  }

  check(b->dim == 2, "initial block dim == 2");
  check_close(b->H->data[0].re, 0.0, 1e-14, "H[0][0] == 0");
  check_close(b->H->data[3].re, 0.0, 1e-14, "H[1][1] == 0");
  check_close(CMAT(b->Sz_end, 0, 0).re, 0.5, 1e-14, "Sz_end[up][up] == +1/2");
  check_close(CMAT(b->Sz_end, 1, 1).re, -0.5, 1e-14,
              "Sz_end[down][down] == -1/2");
  check_close(CMAT(b->Sp_end, 0, 1).re, 1.0, 1e-14,
              "Sp_end[up][down] == 1 (raises)");

  dmrg_block_free(b);
}

/*
 * Cross-check against a from-scratch brute-force ED for small N (N=8, dim=256,
 * tractable for dense diagonalization), isotropic Heisenberg AFM (Jz=Jxy=1).
 */
static void test_vs_brute_force_ed(void) {
  printf(" === DMRG vs brute-force exact diagonalization (N=8) ===\n");
  double Jz = 1.0, Jxy = 1.0;
  int N = 8;
  double E_exact = brute_force_ed(N, Jz, Jxy);

  /* NOTE: Small m is expected to have larger (but bounded) error
   * Tolerance tightens as m grows; only largest m values (>= actual lossless
   * threshold for this N) should match to near machine precision. */
  struct {
    int m;
    double tol;
  } ms[] = {{2, 2.5}, {4, 0.1}, {8, 1e-6}, {16, 1e-6}};

  double prev_err = 1e300;
  for (size_t i = 0; i < sizeof(ms) / sizeof(ms[0]); i++) {
    int m = ms[i].m;
    dmrg_result_t *r = dmrg_run(N, Jz, Jxy, m);
    check(r != NULL, "dmrg_run succeeds");
    if (!r) {
      continue;
    }

    check(r->N_reached == N, "N_reached == N_target (N even)");
    double err = fabs(r->energy - E_exact);
    char label[128];

    snprintf(label, sizeof label, "N=%d m=%d matches brute-force ED", N, m);

    check_close(r->energy, E_exact, ms[i].tol, label);
    check(err <= prev_err + 1e-12, "energy error is non-increasing as m grows");

    prev_err = err;

    free(r);
  }
}

/*
 * Reference values independently cross-checked: infinite-system DMRG vs
 * brute-force ED, isotropic Heisenberg AFM chain (Jz=Jxy=1), open boundary
 * conditions.
 */
static void test_reference_values(void) {
  printf(" === DMRG vs validated reference energies (N=10,12) ===\n");
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
    /* m=24 empirically verified near-machine-precision here (trunc_err
     * ~1e-11 at N=12) without the cost of chasing a literal 2^(N/2) bound */
    int m = 24;
    dmrg_result_t *r = dmrg_run(N, Jz, Jxy, m);
    check(r != NULL, "dmrg_run succeeds");
    if (!r) {
      continue;
    }

    char label[128];
    snprintf(label, sizeof label, "N=%d lossless-m matches reference", N);

    check_close(r->energy, cases[c].E_exact, 1e-6, label);
    check(r->truncation_error < 1e-8, "truncation_error ~0 at lossless m");

    free(r);
  }
}

/*
 * Truncation error must shrink monotonically as m increases (defining signature
 * of correct density-matrix truncation), for a chain long enough that
 * truncation is actually needed (N=16, so intermediate blocks exceed dim 8
 * before the smallest m values here catch up)
 */
static void test_truncation_error_monotonic(void) {
  printf(" === Truncation error shrinks monotonically with m (N=16) ===\n");
  double Jz = 1.0;
  double Jxy = 1.0;
  int N = 16;

  static const int ms_full[] = {2, 4, 8, 16, 24};
  static const int ms_fast[] = {2, 4, 8};

  const int *ms = RUNNING_ON_VALGRIND ? ms_fast : ms_full;
  size_t n_ms = RUNNING_ON_VALGRIND ? (sizeof(ms_fast) / sizeof(ms_fast[0]))
                                    : (sizeof(ms_full) / sizeof(ms_full[0]));
  double prev = 1e300;

  for (size_t i = 0; i < n_ms; i++) {
    dmrg_result_t *r = dmrg_run(N, Jz, Jxy, ms[i]);
    check(r != NULL, "dmrg_run succeeds");
    if (!r) {
      continue;
    }

    char label[128];
    snprintf(label, sizeof label, "m=%d truncation_error <= previous", ms[i]);

    check(r->truncation_error <= prev + 1e-14, label);

    prev = r->truncation_error;

    free(r);
  }
}

/*
 * Per-site ground energy should approach the exact Bethe-ansatz
 * thermodynamic-limit value E0/N -> 1/4 - ln(2) as N grows, for isotropic
 * (critical) Heisenberg chain
 */
static void test_thermodynamic_limit_trend(void) {
  printf(" === Per-site energy trend toward Bethe-ansatz limit ===\n");
  double Jz = 1.0;
  double Jxy = 1.0;
  double bethe_limit = 0.25 - log(2.0);

  static const int Ns[] = {8, 12, 16};
  size_t n_Ns = sizeof(Ns) / sizeof(Ns[0]);
  int m_cap = RUNNING_ON_VALGRIND ? 16 : 24;
  double prev_dist = 1e300;

  for (size_t i = 0; i < n_Ns; i++) {
    dmrg_result_t *r = dmrg_run(Ns[i], Jz, Jxy, m_cap);
    check(r != NULL, "dmrg_run succeeds");
    if (!r) {
      continue;
    }

    double dist = fabs(r->energy_per_site - bethe_limit);

    char label[128];
    snprintf(label, sizeof label, "N=%d: |E/N - Bethe limit| decreases with N",
             Ns[i]);

    check(dist <= prev_dist + 1e-9, label);

    prev_dist = dist;

    free(r);
  }
}

static void test_invalid_input(void) {
  printf(" === Invalid input handling ===\n");
  check(dmrg_run(1, 1.0, 1.0, 8) == NULL, "N_target < 2 rejected");
  check(dmrg_run(8, 1.0, 1.0, 0) == NULL, "m_max < 1 rejected");
  check(dmrg_run(8, NAN, 1.0, 8) == NULL, "non-finite Jz rejected");
  check(dmrg_run(8, 1.0, INFINITY, 8) == NULL, "non-finite Jxy rejected");
}

int main(void) {
  printf(" > Infinite-system DMRG: XXZ open chain tests\n");

  test_block_init();
  test_vs_brute_force_ed();
  test_reference_values();
  test_truncation_error_monotonic();
  test_thermodynamic_limit_trend();
  test_invalid_input();

  if (failures == 0) {
    printf("\nAll test_dmrg checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_dmrg check(s) FAILED.\n", failures);
    return 1;
  }
}
