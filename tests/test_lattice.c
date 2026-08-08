/*
Test: Tight-binding lattice models (1D/2D chains, Anderson localization, SSH
topological edge states).

1. 1D chain (open + periodic): numerical diagonalization must match exact
   analytic dispersion (Bloch's theorem / standing-wave quantization) to near
   machine precision.
2. 2D square lattice (open + periodic): using exact separability of 2D problem
   into independent x/y 1D dispersions.
3. lattice_ipr: deterministic checks on a delta-function state (IPR=1, maximally
   localized) and a uniform superposition (IPR=1/n, maximally delocalized).
4. Anderson localization: disorder_W=0 must reproduce clean chain exactly
   (regression against lattice_1d_chain_analytic); average IPR over several
   disorder realizations must be noticeably larger at strong disorder than at
   weak disorder (qualitative Anderson-localization signature).
5. SSH model: topological phase (t2>t1, open) must show near-zero-energy states
   with large edge weight/IPR; trivial phase (t1>t2) must not.
*/

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/lattice.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

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

static double *diagonalize_sorted(cmatrix_t *H, int n) {
  cmatrix_t *copy = cmatrix_copy(H);
  eigen_t *eig = cmatrix_eigh_complex(copy);
  cmatrix_free(copy);

  double *E = malloc((size_t)n * sizeof *E);
  for (int i = 0; i < n; i++) {
    E[i] = eig->eigenvalues[i];
  }

  eigen_free(eig);

  return E;
}

static void test_1d_chain_vs_analytic(void) {
  printf("test_1d_chain_vs_analytic:\n");

  int n_sites = 10;
  double eps0 = 0.5, t = 1.2;

  for (int bc = 0; bc <= 1; bc++) {
    cmatrix_t *H = lattice_build_1d_chain(n_sites, eps0, t, (lattice_bc_t)bc);
    double *E_num = diagonalize_sorted(H, n_sites);

    double *E_analytic = malloc((size_t)n_sites * sizeof *E_analytic);
    lattice_1d_chain_analytic(n_sites, eps0, t, (lattice_bc_t)bc, E_analytic);

    for (int i = 0; i < n_sites; i++) {
      char label[64];
      snprintf(label, sizeof label, "%s E[%d]",
               bc == LATTICE_PERIODIC ? "periodic" : "open", i);
      check_close(E_num[i], E_analytic[i], 1e-9, label);
    }

    free(E_num);
    free(E_analytic);
    cmatrix_free(H);
  }
}

static void test_2d_square_vs_analytic(void) {
  printf("test_2d_square_vs_analytic:\n");

  int nx = 4, ny = 5;
  double eps0 = 0.3, t = 1.1;
  int N = nx * ny;

  for (int bc = 0; bc <= 1; bc++) {
    cmatrix_t *H = lattice_build_2d_square(nx, ny, eps0, t, (lattice_bc_t)bc);
    double *E_num = diagonalize_sorted(H, N);

    double *E_analytic = malloc((size_t)N * sizeof *E_analytic);
    lattice_2d_square_analytic(nx, ny, eps0, t, (lattice_bc_t)bc, E_analytic);

    double max_err = 0.0;
    for (int i = 0; i < N; i++) {
      double err = fabs(E_num[i] - E_analytic[i]);
      if (err > max_err) {
        max_err = err;
      }
    }

    char label[64];
    snprintf(label, sizeof label,
             "%s: max |E_num - E_analytic| over all %d "
             "states",
             bc == LATTICE_PERIODIC ? "periodic" : "open", N);
    check_close(max_err, 0.0, 1e-8, label);

    free(E_num);
    free(E_analytic);
    cmatrix_free(H);
  }
}

static void test_ipr_deterministic(void) {
  printf("test_ipr_deterministic:\n");

  int n = 8;
  cvector_t *delta = cvector_alloc(n);
  for (int i = 0; i < n; i++) {
    delta->data[i] = c_zero();
  }

  delta->data[3] = c_real(1.0);
  check_close(lattice_ipr(delta), 1.0, 1e-12,
              "delta function: IPR = 1 (maximally localized)");
  cvector_free(delta);

  cvector_t *uniform = cvector_alloc(n);
  for (int i = 0; i < n; i++) {
    uniform->data[i] = c_real(1.0); /* lattice_ipr normalizes internally */
  }
  check_close(lattice_ipr(uniform), 1.0 / n, 1e-12,
              "uniform superposition: IPR = 1/n (maximally delocalized)");
  cvector_free(uniform);

  check_close(lattice_ipr(NULL), 0.0, 0.0, "NULL psi returns 0.0 guard");
}

static void test_anderson_clean_limit(void) {
  printf("test_anderson_clean_limit:\n");

  int n_sites = 12;
  double t = 1.0;

  cmatrix_t *H_anderson = lattice_build_anderson_1d(n_sites, t, 0.0, 42ULL);
  double *E_anderson = diagonalize_sorted(H_anderson, n_sites);

  double *E_clean = malloc((size_t)n_sites * sizeof *E_clean);
  lattice_1d_chain_analytic(n_sites, 0.0, t, LATTICE_OPEN, E_clean);

  double max_err = 0.0;
  for (int i = 0; i < n_sites; i++) {
    double err = fabs(E_anderson[i] - E_clean[i]);
    if (err > max_err) {
      max_err = err;
    }
  }

  check_close(max_err, 0.0, 1e-9,
              "disorder_W=0 Anderson chain matches clean open chain exactly");

  free(E_anderson);
  free(E_clean);
  cmatrix_free(H_anderson);
}

static void test_anderson_localization_trend(void) {
  printf("test_anderson_localization_trend:\n");

  int n_sites = 60;
  double t = 1.0;
  int n_realizations = 15;

  double W_values[2] = {0.2, 6.0}; // weak vs strong disorder
  double avg_ipr[2] = {0.0, 0.0};

  for (int w = 0; w < 2; w++) {
    double sum_ipr = 0.0;

    for (int real = 0; real < n_realizations; real++) {
      cmatrix_t *H = lattice_build_anderson_1d(n_sites, t, W_values[w],
                                               1000ULL + (uint64_t)real);
      cmatrix_t *copy = cmatrix_copy(H);
      eigen_t *eig = cmatrix_eigh_complex(copy);
      cmatrix_free(copy);

      // Using ground state (lowest-energy eigenvector) as localization
      // diagnostic for realization
      cvector_t *psi = cvector_alloc(n_sites);
      for (int i = 0; i < n_sites; i++) {
        psi->data[i] = CMAT(eig->eigenvectors, i, 0);
      }
      sum_ipr += lattice_ipr(psi);

      cvector_free(psi);
      eigen_free(eig);
      cmatrix_free(H);
    }

    avg_ipr[w] = sum_ipr / n_realizations;
  }

  printf("  weak disorder (W=%.1f): avg. ground-state IPR = %.5f\n",
         W_values[0], avg_ipr[0]);
  printf("  strong disorder (W=%.1f): avg. ground-state IPR = %.5f\n",
         W_values[1], avg_ipr[1]);

  check_true(avg_ipr[1] > avg_ipr[0],
             "strong disorder localizes the ground state more than weak "
             "disorder (Anderson localization)");
}

static void test_ssh_topological_edge_states(void) {
  printf("test_ssh_topological_edge_states:\n");

  int n_cells = 15; // 30 sites
  int N = 2 * n_cells;

  // Trivial phase: t1 > t2, no edge states expected
  {
    cmatrix_t *H = lattice_build_ssh(n_cells, 1.0, 0.3, LATTICE_OPEN);
    double *E = diagonalize_sorted(H, N);

    double min_abs_E = 1e9;
    for (int i = 0; i < N; i++) {
      if (fabs(E[i]) < min_abs_E) {
        min_abs_E = fabs(E[i]);
      }
    }

    printf("  trivial (t1=1.0,t2=0.3): closest-to-zero |E| = %.6f\n",
           min_abs_E);
    check_true(min_abs_E > 0.3,
               "trivial phase: no near-zero-energy state (gapped)");

    free(E);
    cmatrix_free(H);
  }

  // Topological phase: t2 > t1, expect two near-zero, edge-localized states
  {
    cmatrix_t *H = lattice_build_ssh(n_cells, 0.3, 1.0, LATTICE_OPEN);
    cmatrix_t *copy = cmatrix_copy(H);
    eigen_t *eig = cmatrix_eigh_complex(copy);
    cmatrix_free(copy);

    // Find eigenvalue closest to zero
    int best = 0;
    double best_abs = fabs(eig->eigenvalues[0]);
    for (int i = 1; i < N; i++) {
      if (fabs(eig->eigenvalues[i]) < best_abs) {
        best_abs = fabs(eig->eigenvalues[i]);
        best = i;
      }
    }

    printf("  topological (t1=0.3,t2=1.0): closest-to-zero E = %.6e\n",
           eig->eigenvalues[best]);
    check_true(best_abs < 1e-3,
               "topological phase: a near-zero-energy state exists");

    cvector_t *psi = cvector_alloc(N);
    for (int i = 0; i < N; i++) {
      psi->data[i] = CMAT(eig->eigenvectors, i, best);
    }

    double edge_weight = c_abs2(psi->data[0]) + c_abs2(psi->data[1]) +
                         c_abs2(psi->data[N - 1]) + c_abs2(psi->data[N - 2]);
    double ipr = lattice_ipr(psi);

    printf("  edge weight (first+last 2 sites) = %.4f   IPR = %.4f\n",
           edge_weight, ipr);
    check_true(edge_weight > 0.5,
               "topological zero mode is concentrated at the chain edges");
    check_true(ipr > 0.1, "topological zero mode is significantly more "
                          "localized than a generic extended state (~1/N)");

    cvector_free(psi);
    eigen_free(eig);
    cmatrix_free(H);
  }
}

static void test_landau_zero_field_matches_plain_square(void) {
  printf("test_landau_zero_field_matches_plain_square:\n");

  int nx = 6, ny = 6, N = nx * ny;
  double eps0 = 0.3, t = 1.1;

  // NOTE: \alpha=0 must reduce lattice_build_2d_square_magnetic exactly to
  // lattice_build_2d_square (all Peierls phases = \exp(0) = 1). Using
  // LATTICE_OPEN on both sides for a like-for-like comparison: the magnetic
  // builder is always open in x (required by the Landau gauge choice), so
  // comparing against a plain square lattice that is periodic in x would differ
  // by x-wraparound bonds alone
  cmatrix_t *H_mag =
      lattice_build_2d_square_magnetic(nx, ny, eps0, t, 0.0, LATTICE_OPEN);
  cmatrix_t *H_plain = lattice_build_2d_square(nx, ny, eps0, t, LATTICE_OPEN);

  double max_entry_err = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      complex_t a = CMAT(H_mag, i, j);
      complex_t b = CMAT(H_plain, i, j);
      double err = c_abs(c_sub(a, b));
      if (err > max_entry_err) {
        max_entry_err = err;
      }
    }
  }
  check_close(max_entry_err, 0.0, 1e-12,
              "\\alpha=0 magnetic-lattice matrix == plain square-lattice "
              "matrix, entrywise");

  double *E_num = diagonalize_sorted(H_mag, N);
  double *E_analytic = malloc((size_t)N * sizeof *E_analytic);
  lattice_2d_square_analytic(nx, ny, eps0, t, LATTICE_OPEN, E_analytic);

  double max_eig_err = 0.0;
  for (int i = 0; i < N; i++) {
    double err = fabs(E_num[i] - E_analytic[i]);
    if (err > max_eig_err) {
      max_eig_err = err;
    }
  }
  check_close(max_eig_err, 0.0, 1e-8,
              "\\alpha=0 eigenvalues match plain-square analytic dispersion");

  free(E_num);
  free(E_analytic);
  cmatrix_free(H_mag);
  cmatrix_free(H_plain);
}

static void test_landau_hermiticity(void) {
  printf("test_landau_hermiticity:\n");

  int nx = 10, ny = 10, N = nx * ny;
  cmatrix_t *H = lattice_build_2d_square_magnetic(nx, ny, 0.0, 1.0, 0.037,
                                                  LATTICE_PERIODIC);

  double max_err = 0.0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      complex_t err = c_sub(CMAT(H, i, j), c_conj(CMAT(H, j, i)));
      double e = c_abs(err);
      if (e > max_err) {
        max_err = e;
      }
    }
  }
  check_close(max_err, 0.0, 1e-12, "H == H^dagger entrywise (alpha != 0)");

  cmatrix_free(H);
}

static void test_landau_continuum_limit(void) {
  printf("test_landau_continuum_limit:\n");

  // NOTE: Weak-field limit: ground Landau level should sit close to
  // continuum-limit prediction on a lattice much larger than magnetic length ~
  // 1/\sqrt(2 * \pi * \alpha), away from open-boundary edge-state
  // contamination.
  double t = 1.0, alpha = 0.05;
  int nx = 10, ny = 10, N = nx * ny;

  cmatrix_t *H =
      lattice_build_2d_square_magnetic(nx, ny, 0.0, t, alpha, LATTICE_PERIODIC);
  double *E = diagonalize_sorted(H, N);

  double predicted_n0 = lattice_landau_level_energy(0, 0.0, t, alpha);
  double rel_err = fabs(E[0] - predicted_n0) / fabs(predicted_n0);

  printf("  n=0: numeric=%.6f predicted=%.6f rel_err=%.4f%%\n", E[0],
         predicted_n0, rel_err * 100.0);
  check_true(rel_err < 0.02,
             "n=0 Landau level within 2% of continuum-limit prediction");

  free(E);
  cmatrix_free(H);
}

int main(void) {
  test_1d_chain_vs_analytic();
  test_2d_square_vs_analytic();
  test_ipr_deterministic();
  test_anderson_clean_limit();
  test_anderson_localization_trend();
  test_ssh_topological_edge_states();
  test_landau_zero_field_matches_plain_square();
  test_landau_hermiticity();
  test_landau_continuum_limit();

  if (failures == 0) {
    printf("\nAll test_lattice checks passed.\n");
    return 0;
  } else {
    printf("\n%d check(s) FAILED.\n", failures);
    return 1;
  }
}
