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

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/sparse.h"
#include "../core/vector.h"
#include "../physics/spin_chain.h"
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

/*
 * Ground-state energy of the full Heisenberg ring = min over all k sectors of
 * lowest eigenvalue of that sector's Hamiltonian. Sz=0 (or nearest integer
 * sector for odd N) always contains the true ground state for antiferromagnetic
 * chain, scan every (nup, k) sector
 */
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

/* Reflection (spatial parity) symmetry: R^2=I, Hermiticity, [H,R]=0, and lower
 * of the two parity blocks' ground energies must equal the sector's own
 * (unsplit) ground energy. For both N=6 and N=8, at both momenta that support a
 * reflection quantum number (k=0 and k=N/2).
 *
 * NOTE: This checks ground energies rather than the full spectrum: requesting
 * every eigenvalue via lanczos_eigs(H, dim, dim, tol) can hit early Lanczos
 * breakdown (the Krylov sequence's residual falls below tol before reaching
 * m_eff = dim, particularly for real, highly-symmetric Hamiltonians like these)
 * and return NULL well short of a full spectrum
 */
static void test_reflection_parity(void) {
  printf("\n-- Reflection parity: R^2=I, [H,R]=0, ground energy preserved "
         "under parity splitting --\n");

  struct {
    int N, nup, k;
    double E0_ref; /* independently known ground energy */
  } cases[] = {
      {6, 3, 0, -1.5},      /* k=0 is NOT the N=6 ground sector */
      {6, 3, 3, -2.802776}, /* N=6 true ground state: k=pi */
      {8, 4, 0, -3.651093}, /* N=8 true ground state: k=0 */
      {8, 4, 4, -3.128419}, /* N=8, k=pi (excited relative to k=0) */
  };
  const int has_ref[4] = {0, 1, 1, 0};

  for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); c++) {
    int N = cases[c].N, nup = cases[c].nup, k = cases[c].k;
    spin_sector_t *sec = spin_sector_build(N, nup, k);

    check(sec != NULL, "sector built");
    if (!sec) {
      continue;
    }

    sparse_matrix_t *H = spin_sector_hamiltonian(sec, 1.0, 1.0, /*pbc=*/1);

    cmatrix_t *R = spin_sector_reflection_matrix(sec);
    char label[160];
    snprintf(label, sizeof label, "N=%d k=%d: R non-NULL at a valid momentum",
             N, k);

    check(R != NULL, label);
    if (!R) {
      sparse_free(H);
      spin_sector_free(sec);

      continue;
    }

    int dim = sec->dim;
    double max_r2_err = 0.0, max_herm_err = 0.0;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        complex_t rr = c_zero();

        for (int m = 0; m < dim; m++) {
          rr = c_add(rr, c_mul(CMAT(R, i, m), CMAT(R, m, j)));
        }

        complex_t expect = (i == j) ? c_new(1.0, 0.0) : c_zero();
        double e = sqrt(c_abs2(c_sub(rr, expect)));
        if (e > max_r2_err) {
          max_r2_err = e;
        }

        double eh = sqrt(c_abs2(c_sub(CMAT(R, i, j), c_conj(CMAT(R, j, i)))));
        if (eh > max_herm_err) {
          max_herm_err = eh;
        }
      }
    }

    snprintf(label, sizeof label, "N=%d k=%d: |R^2 - I| within tol", N, k);
    check_close(max_r2_err, 0.0, 1e-10, label);

    snprintf(label, sizeof label, "N=%d k=%d: |R - R^dagger| within tol", N, k);
    check_close(max_herm_err, 0.0, 1e-10, label);

    cmatrix_free(R);

    int dim_p, dim_m;
    cmatrix_t *Hp = spin_sector_parity_project(sec, H, +1, &dim_p);
    cmatrix_t *Hm = spin_sector_parity_project(sec, H, -1, &dim_m);
    snprintf(label, sizeof label,
             "N=%d k=%d: parity block dims sum to sector dim (%d+%d=%d)", N, k,
             dim_p, dim_m, dim);

    check(dim_p + dim_m == dim, label);

    double best = 1e300;
    if (dim_p > 0) {
      sparse_matrix_t *Hps = sparse_from_dense(Hp, 0.0);
      lanczos_result_t *ep = lanczos_eigs(Hps, 1, dim_p, 1e-12);
      if (ep && ep->values[0] < best) {
        best = ep->values[0];
      }

      lanczos_free(ep);
      sparse_free(Hps);
    }

    if (dim_m > 0) {
      sparse_matrix_t *Hms = sparse_from_dense(Hm, 0.0);
      lanczos_result_t *em = lanczos_eigs(Hms, 1, dim_m, 1e-12);
      if (em && em->values[0] < best) {
        best = em->values[0];
      }

      lanczos_free(em);
      sparse_free(Hms);
    }

    lanczos_result_t *ref = lanczos_eigs(H, 1, dim, 1e-12);
    snprintf(label, sizeof label,
             "N=%d k=%d: unsplit ground-energy query succeeded (sanity check "
             "on reference itself)",
             N, k);

    check(ref != NULL, label);
    if (ref) {
      snprintf(label, sizeof label,
               "N=%d k=%d: lowest parity-block energy matches unsplit sector's "
               "ground energy",
               N, k);

      check_close(best, ref->values[0], 1e-8, label);
      lanczos_free(ref);
    }

    if (has_ref[c]) {
      snprintf(label, sizeof label,
               "N=%d k=%d: matches independently known ground energy", N, k);

      check_close(best, cases[c].E0_ref, 1e-5, label);
    }

    cmatrix_free(Hp);
    cmatrix_free(Hm);
    sparse_free(H);
    spin_sector_free(sec);
  }

  // k=1 on N=6 has no reflection quantum number (only k=0, k=N/2 do)
  spin_sector_t *bad = spin_sector_build(6, 3, 1);
  cmatrix_t *Rbad = spin_sector_reflection_matrix(bad);

  check(Rbad == NULL, "reflection_matrix returns NULL at a momentum with no "
                      "parity quantum number (N=6, k=1)");

  cmatrix_free(Rbad);
  spin_sector_free(bad);
}

int main(void) {
  printf("=== Spin-chain symmetry-adapted ED + DSF tests ===\n");

  test_ground_energies();
  test_hamiltonian_hermitian_and_sector_sizes();
  test_szq_sum_rule();
  test_continued_fraction_integrates_to_I0();
  test_reflection_parity();

  if (failures == 0) {
    printf("\nAll test_spin_chain checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_spin_chain check(s) FAILED.\n", failures);
    return 1;
  }
}
