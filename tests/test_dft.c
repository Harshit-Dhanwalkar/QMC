/*
 * Test: closed-shell Kohn-Sham LDA DFT (Slater exchange + PZ81 correlation) for
 * s-orbitals-only atoms/ions.
 *
 * Same grid-resolution as test_hartree_fock.c:
 * finite-difference/dense-diagonalization solve on a finite radial grid, not a
 * converged calculation, so tolerances are meaningful but grid-tolerant, not
 * tight to several decimal places.
 *
 * Cross-checked against References: RHF limit, Szabo & Ostlund, Pekeris (1959),
 * Hartree and arXiv:2202.00647, Table II
 */

#include "../core/utils.h"
#include "../core/vector.h"
#include "../physics/dft.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef RUNNING_ON_VALGRIND
#define RUNNING_ON_VALGRIND 0
#endif

#define HELIUM_EXACT_NONREL -2.903724
#define HELIUM_HF -2.861680
#define HELIUM_VWN_LDA_REFERENCE -2.834836

static int failures = 0;

static int check_range(double got, double lo, double hi, const char *label) {
  printf("  %s: got=%.6f  expected range=[%.6f, %.6f]\n", label, got, lo, hi);

  return (got < lo || got > hi);
}

static void check_true(int cond, const char *label) {
  printf("  %s: %s\n", label, cond ? "ok" : "FAIL");
  if (!cond) {
    failures++;
  }
}

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

/* Both V_x and V_c_pz81 must equal thermodynamic relation:
 * V_xc = d(n * eps_xc) / dn, checked via central finite difference
 */
static void test_xc_potential_matches_finite_difference(void) {
  printf("test_xc_potential_matches_finite_difference:\n");

  const double densities[] = {0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0};
  for (int i = 0; i < 7; i++) {
    double n = densities[i];
    double h = n * 1e-6;

    double fd_x = ((n + h) * lda_exchange_energy_density(n + h) -
                   (n - h) * lda_exchange_energy_density(n - h)) /
                  (2.0 * h);
    double fd_c = ((n + h) * lda_correlation_energy_density_pz81(n + h) -
                   (n - h) * lda_correlation_energy_density_pz81(n - h)) /
                  (2.0 * h);

    char label_x[64], label_c[64];
    snprintf(label_x, sizeof(label_x), "V_x(%.3f) matches d(n * eps_x)/dn", n);
    snprintf(label_c, sizeof(label_c), "V_c(%.3f) matches d(n * eps_c)/dn", n);

    check_close(lda_exchange_potential(n), fd_x, 1e-6, label_x);
    check_close(lda_correlation_potential_pz81(n), fd_c, 1e-6, label_c);
  }
}

static void test_xc_sign_and_limits(void) {
  printf("test_xc_sign_and_limits:\n");

  check_true(lda_exchange_potential(1.0) < 0.0,
             "exchange potential is attractive (negative)");
  check_true(lda_correlation_potential_pz81(1.0) < 0.0,
             "correlation potential is attractive (negative)");
  check_true(lda_exchange_potential(0.0) == 0.0,
             "exchange potential vanishes at zero density");
  check_true(lda_correlation_potential_pz81(0.0) == 0.0,
             "correlation potential vanishes at zero density");

  // |V_x| must increase with density (n^{1/3} is monotonic increasing)
  check_true(fabs(lda_exchange_potential(2.0)) >
                 fabs(lda_exchange_potential(1.0)),
             "|V_x| increases with density");
}

static void test_helium_ks_lda(void) {
  printf("test_helium_ks_lda:\n");

  int N = RUNNING_ON_VALGRIND ? 40 : 160;
  double r_min = 1e-4;
  double r_max = RUNNING_ON_VALGRIND ? 6.0 : 12.0;
  double *r = linspace(r_min, r_max, N);
  int max_iter = RUNNING_ON_VALGRIND ? 50 : 200;

  dft_result_t *res =
      dft_lda_atom_s_orbitals(r, N, 2.0, 1, 0.3, 1e-9, max_iter);
  check_true(res != NULL, "dft_lda_atom_s_orbitals allocates for He");
  if (!res) {
    free(r);
    return;
  }

  printf("  He KS-LDA: E_total=%.6f  eps_1s=%.6f  converged=%d  iters=%d\n",
         res->total_energy, res->orbital_energies[0], res->converged,
         res->iterations);
  printf("  Reference points: exact=%.6f  HF=%.6f  VWN-LDA=%.6f\n",
         HELIUM_EXACT_NONREL, HELIUM_HF, HELIUM_VWN_LDA_REFERENCE);

  check_true(res->converged, "He SCF converges within max_iter");

  // Range relaxed under Valgrind (coarser grid → less accurate)
  double lo = RUNNING_ON_VALGRIND ? -2.9 : -2.87;
  double hi = RUNNING_ON_VALGRIND ? -2.5 : -2.70;
  check_true(!check_range(res->total_energy, lo, hi,
                          "He total_energy in physically sane range"),
             "He total_energy in physically sane range");

  double dr = r[1] - r[0];
  double norm_sq = 0.0;
  for (int i = 0; i < N; i++) {
    double u = res->orbitals[0]->data[i].re;
    norm_sq += u * u * dr;
  }

  check_close(norm_sq, 1.0, 1e-6,
              "occupied orbital normalized: \\int u^2 dr = 1");

  check_true(res->E_hartree > 0.0,
             "E_hartree (classical e-e repulsion) is positive");
  check_true(res->E_xc < 0.0, "E_xc is negative (exchange dominates)");

  dft_result_free(res);
  free(r);
}

static void test_invalid_input_rejected(void) {
  printf("test_invalid_input_rejected:\n");

  double r[20];
  for (int i = 0; i < 20; i++) {
    r[i] = 1e-4 + i * 0.1;
  }

  check_true(dft_lda_atom_s_orbitals(NULL, 20, 2.0, 1, 0.3, 1e-9, 50) == NULL,
             "NULL grid rejected");
  check_true(dft_lda_atom_s_orbitals(r, 5, 2.0, 1, 0.3, 1e-9, 50) == NULL,
             "N < 10 rejected");
  check_true(dft_lda_atom_s_orbitals(r, 20, -1.0, 1, 0.3, 1e-9, 50) == NULL,
             "Z <= 0 rejected");
  check_true(dft_lda_atom_s_orbitals(r, 20, 2.0, 0, 0.3, 1e-9, 50) == NULL,
             "n_orbitals < 1 rejected");
  check_true(dft_lda_atom_s_orbitals(r, 20, 2.0, 1, 1.5, 1e-9, 50) == NULL,
             "mix out of (0,1] rejected");
}

int main(void) {
  test_xc_potential_matches_finite_difference();
  test_xc_sign_and_limits();
  test_invalid_input_rejected();
  test_helium_ks_lda();

  if (failures == 0) {
    printf("\nAll test_dft checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_dft check(s) FAILED.\n", failures);
    return 1;
  }
}
