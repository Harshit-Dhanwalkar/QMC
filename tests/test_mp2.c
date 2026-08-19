/*
 * Test: MP2 (s-orbitals-only restricted) second-order perturbation theory.
 *
 * 1. Deterministic fixture: mp2_correlation_energy on a built hf_result_t
 *    with synthetic (non-HF) occupied/virtual orbitals and orbital energies,
 * 2. Physics: MP2 on top of converged helium HF must be negative, monotonically
 *    increase (in magnitude) as more virtuals are included, and move total
 *    energy toward the exact ground state.
 *     NOTE: Given the s-only restriction , this recovers only part of true
 *    correlation energy (angular correlation from p/d virtuals is missing), so
 *    test currently checks direction, not exact recovery.
 *     FIX: Now QMC support p and d orbitals
 * 3. Invalid-input handling.
 */

#include "../core/complex.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../physics/hartree_fock.h"
#include "../physics/mp2.h"
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
  printf("  %s: got=%.10f expected=%.10f err=%.2e\n", label, got, expected,
         err);
  if (err > tol) {
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

static double *make_radial(const double *r, int N, double dr, double a,
                           double b) {
  double *f = malloc((size_t)N * sizeof *f);
  double norm_sq = 0.0;
  for (int i = 0; i < N; i++) {
    f[i] = r[i] * exp(-a * r[i]) + 0.3 * r[i] * r[i] * exp(-b * r[i]);
    norm_sq += f[i] * f[i] * dr;
  }

  double norm = sqrt(norm_sq);
  for (int i = 0; i < N; i++) {
    f[i] /= norm;
  }

  return f;
}

// Builds a synthetic hf_result_t with 2 occupied and 2 virtual orbitals
static hf_result_t *build_synthetic_hf(const double *r, int N) {
  hf_result_t *hf = malloc(sizeof *hf);

  hf->n_orbitals = 2;
  hf->N = N;
  hf->Z = 0.0;
  hf->total_energy = 0.0;
  hf->iterations = 0;
  hf->converged = 1;

  double dr = r[1] - r[0];

  hf->orbital_energies = malloc(2 * sizeof(double));
  hf->orbital_energies[0] = -1.2;
  hf->orbital_energies[1] = -0.9;

  hf->orbitals = malloc(2 * sizeof(cvector_t *));
  double *o0 = make_radial(r, N, dr, 1.0, 0.8);
  double *o1 = make_radial(r, N, dr, 0.6, 1.4);

  hf->orbitals[0] = cvector_alloc(N);
  hf->orbitals[1] = cvector_alloc(N);

  for (int i = 0; i < N; i++) {
    hf->orbitals[0]->data[i] = c_real(o0[i]);
    hf->orbitals[1]->data[i] = c_real(o1[i]);
  }

  free(o0);
  free(o1);

  hf->n_virtual = 2;
  hf->virtual_energies = malloc(2 * sizeof(double));
  hf->virtual_energies[0] = 0.3;
  hf->virtual_energies[1] = 0.5;

  hf->virtual_orbitals = malloc(2 * sizeof(cvector_t *));

  double *v0 = make_radial(r, N, dr, 1.5, 1.2);
  double *v1 = make_radial(r, N, dr, 2.0, 0.6);

  hf->virtual_orbitals[0] = cvector_alloc(N);
  hf->virtual_orbitals[1] = cvector_alloc(N);

  for (int i = 0; i < N; i++) {
    hf->virtual_orbitals[0]->data[i] = c_real(v0[i]);
    hf->virtual_orbitals[1]->data[i] = c_real(v1[i]);
  }

  free(v0);
  free(v1);

  return hf;
}

// Deterministic fixture must use same grid as the reference value.
static void test_synthetic_fixture(void) {
  printf("test_synthetic_fixture:\n");

  int N = 160;
  double *r = linspace(1e-4, 12.0, N);
  hf_result_t *hf = build_synthetic_hf(r, N);

  mp2_result_t res = mp2_correlation_energy(hf, r, N, 2);

  check_close(res.e_mp2, -0.6561117086134629, 1e-6,
              "E_MP2 matches reference implementation");
  check_true(res.n_occ == 2 && res.n_virt == 2,
             "n_occ/n_virt reported correctly");

  hf_result_free(hf);
  free(r);
}

static void test_helium_mp2_physical_sanity(void) {
  printf("test_helium_mp2_physical_sanity:\n");

  int N = RUNNING_ON_VALGRIND ? 40 : 160;
  double r_min = 1e-4;
  double r_max = RUNNING_ON_VALGRIND ? 6.0 : 12.0;
  double *r = linspace(r_min, r_max, N);
  int max_iter = RUNNING_ON_VALGRIND ? 50 : 200;

  hf_result_t *hf =
      hartree_fock_atom_s_orbitals(r, N, 2.0, 1, 0.4, 1e-8, max_iter);
  check_true(hf != NULL && hf->converged, "helium HF reference converges");

  if (hf) {
    double E_exact = -2.9037;
    double prev_mp2_magnitude = 0.0;
    const int nv_list[4] = {2, 5, 10, 20};
    int nv_count = RUNNING_ON_VALGRIND ? 2 : 4;

    for (int t = 0; t < nv_count; t++) {
      mp2_result_t res = mp2_correlation_energy(hf, r, N, nv_list[t]);

      printf("  n_virtual=%2d: E_MP2=%.8f  E_HF+MP2=%.6f\n", nv_list[t],
             res.e_mp2, res.e_total);

      char label[64];
      snprintf(label, sizeof label, "n_virtual=%d: E_MP2 is negative",
               nv_list[t]);
      check_true(res.e_mp2 < 0.0, label);

      snprintf(label, sizeof label,
               "n_virtual=%d: |E_MP2| grows monotonically with more virtuals",
               nv_list[t]);
      check_true(fabs(res.e_mp2) >= prev_mp2_magnitude - 1e-12, label);
      prev_mp2_magnitude = fabs(res.e_mp2);

      snprintf(label, sizeof label,
               "n_virtual=%d: E_HF+MP2 closer to exact than E_HF alone",
               nv_list[t]);
      check_true(fabs(res.e_total - E_exact) < fabs(hf->total_energy - E_exact),
                 label);
    }

    hf_result_free(hf);
  }

  free(r);
}

static void test_invalid_input(void) {
  printf("test_invalid_input:\n");

  int N = RUNNING_ON_VALGRIND ? 40 : 160;
  double *r = linspace(1e-4, RUNNING_ON_VALGRIND ? 6.0 : 12.0, N);
  hf_result_t *hf = build_synthetic_hf(r, N);

  mp2_result_t r1 = mp2_correlation_energy(NULL, r, N, 1);
  check_true(r1.e_mp2 == 0.0 && r1.n_occ == 0, "NULL hf rejected");

  mp2_result_t r2 = mp2_correlation_energy(hf, NULL, N, 1);
  check_true(r2.e_mp2 == 0.0 && r2.n_occ == 0, "NULL r rejected");

  mp2_result_t r3 = mp2_correlation_energy(hf, r, N, 0);
  check_true(r3.e_mp2 == 0.0 && r3.n_occ == 0, "n_virtual=0 rejected");

  mp2_result_t r4 = mp2_correlation_energy(hf, r, N, 100);
  check_true(r4.e_mp2 == 0.0 && r4.n_occ == 0,
             "n_virtual > hf->n_virtual rejected");

  hf_result_free(hf);
  free(r);
}

int main(void) {
  test_synthetic_fixture();
  test_invalid_input();
  test_helium_mp2_physical_sanity();

  if (failures == 0) {
    printf("\nAll test_mp2 checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
