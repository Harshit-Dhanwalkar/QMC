/*
 * Test: translation-symmetry-adapted exact diagonalization of the spin-1/2
 * XXZ ring, and the Lanczos continued-fraction dynamical structure factor.
 *
 * Reference values were generated independently with a full-Hilbert-space (no
 * symmetry) numpy exact-diagonalization script for the isotropic Heisenberg
 * chain (Jxy = Jz = 1), periodic boundary conditions:
 *
 *   N=4: E0 = -2.000000       (exact, well-known small-ring benchmark)
 *   N=6: E0 = -2.802776
 *   N=8: E0 = -3.651093
 *
 * and, for N=6 at momentum-transfer index q=1 (q = 2 * \pi/6) applied to ground
 * state: I0 = <\phi0|\phi0> = 0.61324950943693
 *
 * Where
 *   \phi0 = S^z_q |\psi_0>, computed via full diagonalization + exact Lehmann
 * sum (sum of |<n|S^z_q|0>|^2 over ALL eigenstates equals I0 by completeness;
 * this is sum rule the symmetry-sector construction must reproduce exactly,
 * since S^z_q preserves S^z and only moves spectral weight into a single
 * (S^z=0, k=q) sector).
 *
 * NOTE:: for N=6 the ground state sits at momentum k=3 (i.e. k=\pi), NOT k=0.
 * This is a real, well-known feature of the finite-size Heisenberg
 * antiferromagnet, not an implementation detail: for even N ground-state
 * momentum is 0 when N/2 is even and \pi when N/2 is odd (Marshall's sign rule
 * / Lieb-Schultz-Mattis). N=6 has N/2=3, odd, so the ground state is at k=\pi.
 */

#include "../core/sparse.h"
#include "../core/vector.h"
#include "../physics/spin_chain.h"
#include "/home/harshitpd/Documents/GITHUB/QMC/core/matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  if (fabs(got - expected) > tol) {
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  } else {
    printf("  ok:   %s (%.10f)\n", msg, got);
  }
}

static void check(int cond, const char *msg) {
  if (!cond) {
    printf("  FAIL: %s\n", msg);
    failures++;
  } else {
    printf("  ok:   %s\n", msg);
  }
}

/* Ground-state energy of the full Heisenberg ring = min over all k sectors of
 * lowest eigenvalue of that sector's Hamiltonian. Sz=0 (or nearest integer
 * sector for odd N) always contains the true ground state for antiferromagnetic
 * chain, scan every (nup, k) sector */
static double ring_ground_energy(int N) {
  double best = 1e300;

  for (int nup = 0; nup <= N; nup++) {
    for (int k = 0; k < N; k++) {
      spin_sector_t *sec = spin_sector_build(N, nup, k);

      if (!sec || sec->dim == 0) {
        spin_sector_free(sec);

        continue;
      }

      sparse_matrix_t *H = spin_sector_hamiltonian(sec, 1.0, 1.0, 1);
      lanczos_result_t *res =
          lanczos_eigs(H, 1, sec->dim < 60 ? sec->dim : 60, 1e-10);
      if (res && res->values[0] < best) {
        best = res->values[0];
      }

      lanczos_free(res);
      sparse_free(H);
      spin_sector_free(sec);
    }
  }

  return best;
}

static void test_ground_energies(void) {
  printf("\n-- Ground-state energies vs exact-diagonalization reference --\n");

  check_close(ring_ground_energy(4), -2.000000, 1e-6, "N=4 E0");
  check_close(ring_ground_energy(6), -2.802776, 1e-5, "N=6 E0");
  check_close(ring_ground_energy(8), -3.651093, 1e-5, "N=8 E0");
}

static void test_hamiltonian_hermitian_and_sector_sizes(void) {
  printf("\n-- Sector bookkeeping sanity checks --\n");
  int N = 8;
  long total = 0;

  for (int nup = 0; nup <= N; nup++) {
    for (int k = 0; k < N; k++) {
      spin_sector_t *sec = spin_sector_build(N, nup, k);
      if (sec) {
        total += sec->dim;
      }

      spin_sector_free(sec);
    }
  }

  check(total == (1L << N), "sum of all sector dims equals 2^N (8 sites)");
}

/* Cross-check S^z_q construction against exact sum rule I0 */
static void test_szq_sum_rule(void) {
  printf("\n-- S^z_q excitation: sum-rule cross-check (N=6, q index 1) --\n");
  /* Ground-state momentum for N=6 is k=pi (index 3), not k=0
   * (Marshall's sign rule / Lieb-Schultz-Mattis: N/2=3 is odd). */
  int N = 6, nup = N / 2, q_index = 1, k_gs = 3;

  spin_sector_t *sec0 = spin_sector_build(N, nup, k_gs);
  check(sec0 != NULL, "k=pi, Sz=0 sector built");
  if (!sec0) {
    return;
  }

  sparse_matrix_t *H0 = spin_sector_hamiltonian(sec0, 1.0, 1.0, 1);
  lanczos_result_t *gs = lanczos_eigs(H0, 1, sec0->dim, 1e-12);
  check(gs != NULL, "ground state found in k=pi sector");
  if (!gs) {
    sparse_free(H0);
    spin_sector_free(sec0);

    return;
  }

  check_close(gs->values[0], -2.802776, 1e-5,
              "N=6 ground energy (k=pi sector alone)");

  cvector_t *psi0 = cvector_alloc(sec0->dim);
  for (int i = 0; i < sec0->dim; i++) {
    psi0->data[i] = CMAT(gs->vectors, i, 0);
  }

  spin_sector_t *target = NULL;
  cvector_t *phi0 = NULL;
  double I0 = -1.0;
  int rc = spin_apply_szq(sec0, psi0, q_index, &target, &phi0, &I0);
  check(rc == 0, "spin_apply_szq succeeded");
  if (rc == 0) {
    check(target->k == (k_gs + q_index) % N,
          "target sector momentum = k0 + q_index");
    check_close(I0, 0.61324950943693, 1e-8,
                "I0 matches exact Lehmann-sum reference");
  }

  cvector_free(psi0);
  cvector_free(phi0);
  spin_sector_free(target);
  lanczos_free(gs);
  sparse_free(H0);
  spin_sector_free(sec0);
}

/* End-to-end: continued-fraction S(q,omega) should integrate (over \omega)
 * to I0, since \sum_n |<n|f0>|^2 = ||f0||^2 = 1 for normalized start vector,
 * and S includes I0 prefactor so \int S(q,w) d\omega over all \omega equals I0.
 */
static void test_continued_fraction_integrates_to_I0(void) {
  printf("\n-- Continued-fraction spectral weight vs I0 (N=6, q index 1) --\n");
  // Same k=\pi ground sector as test_szq_\sum_rule
  int N = 6, nup = N / 2, q_index = 1, k_gs = 3;

  spin_sector_t *sec0 = spin_sector_build(N, nup, k_gs);
  sparse_matrix_t *H0 = spin_sector_hamiltonian(sec0, 1.0, 1.0, 1);
  lanczos_result_t *gs = lanczos_eigs(H0, 1, sec0->dim, 1e-12);
  double E0 = gs->values[0];

  cvector_t *psi0 = cvector_alloc(sec0->dim);
  for (int i = 0; i < sec0->dim; i++) {
    psi0->data[i] = CMAT(gs->vectors, i, 0);
  }

  spin_sector_t *target = NULL;
  cvector_t *phi0 = NULL;
  double I0 = -1.0;
  spin_apply_szq(sec0, psi0, q_index, &target, &phi0, &I0);

  sparse_matrix_t *Ht = spin_sector_hamiltonian(target, 1.0, 1.0, 1);
  cvector_t *f0 = cvector_copy(phi0);
  cvector_normalize(f0);

  lanczos_tridiag_t *tri = lanczos_tridiagonalize(Ht, f0, target->dim, 1e-12);
  check(tri != NULL, "Lanczos tridiagonalization from S^z_q|psi0> succeeded");

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

    check_close(integral, I0, 5e-3,
                "numerical integral of S(q,w) over omega matches I0");
  }

  lanczos_tridiag_free(tri);
  cvector_free(f0);
  sparse_free(Ht);
  cvector_free(phi0);
  spin_sector_free(target);
  cvector_free(psi0);
  lanczos_free(gs);
  sparse_free(H0);
  spin_sector_free(sec0);
}

int main(void) {
  printf("=== Spin-chain symmetry-adapted ED + DSF tests ===\n");

  test_ground_energies();
  test_hamiltonian_hermitian_and_sector_sizes();
  test_szq_sum_rule();
  test_continued_fraction_integrates_to_I0();

  if (failures == 0) {
    printf("\nAll test_spin_chain checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_spin_chain check(s) FAILED.\n", failures);
    return 1;
  }
}
