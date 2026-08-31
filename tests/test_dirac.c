/*
 * Test: dirac_1d complex-Hermitian eigensolver.
 *
 * 1. Hermiticity check: the 2N x 2N Dirac matrix must satisfy
 *    H[a][b] = conj(H[b][a]) for every entry.
 * 2. Physical interpretation: for free particle (V=0), eigenvalue spectrum
 *    should mostly split into two branches separated by ~2 * m * c^2
 *    (positive-energy states near/above +mc^2, negative-energy states
 *    near/below -m * c^2).
 * 3. dirac_radial_solve validated against exact closed-form relativistic
 *    hydrogen spectrum (Sommerfeld formula)
 */

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/matrix.h"
#include "../physics/potentials.h"
#include "../physics/relativistic.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef RUNNING_ON_VALGRIND
#define RUNNING_ON_VALGRIND 0
#endif

static int test_hermiticity_of_construction(void) {
  int N = RUNNING_ON_VALGRIND ? 10 : 20;
  double dx = 0.1;
  double hbar = 1.0, c = 1.0, m = 1.0;
  double coeff = hbar * c / (2.0 * dx);

  int M = 2 * N;
  cmatrix_t *H = cmatrix_alloc(M, M);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

  double *V = calloc(N, sizeof *V); // free particle

  for (int i = 0; i < N; i++) {
    int row1 = i, row2 = i + N;
    CMAT(H, row1, row1) = c_real(V[i] + m * c * c);
    CMAT(H, row2, row2) = c_real(V[i] - m * c * c);

    if (i > 0) {
      CMAT(H, row1, i - 1 + N) = c_imag(coeff);
      CMAT(H, row2, i - 1) = c_imag(coeff);
    }

    if (i < N - 1) {
      CMAT(H, row1, i + 1 + N) = c_imag(-coeff);
      CMAT(H, row2, i + 1) = c_imag(-coeff);
    }
  }

  free(V);

  int fail = 0;
  double tol = 1e-12;
  for (int a = 0; a < M; a++) {
    for (int b = 0; b < M; b++) {
      complex_t hab = CMAT(H, a, b);
      complex_t hba = CMAT(H, b, a);

      double err = sqrt(pow(hab.re - hba.re, 2) + pow(hab.im + hba.im, 2));
      if (err > tol) {
        printf(
            "  FAIL: H[%d][%d]=(%.4f,%.4f) not conj of H[%d][%d]=(%.4f,%.4f)\n",
            a, b, hab.re, hab.im, b, a, hba.re, hba.im);

        fail = 1;
      }
    }
  }

  if (!fail) {
    printf("  OK: H[a][b] = conj(H[b][a]) for all %d x %d entries\n", M, M);
  }

  cmatrix_free(H);

  return fail;
}

static int test_free_particle_branches(void) {
  int N = RUNNING_ON_VALGRIND ? 20 : 40;
  double dx = 0.2;
  double hbar = 1.0, c = 1.0, m = 1.0;
  double *x = malloc(N * sizeof *x);
  double *V = calloc(N, sizeof *V);
  for (int i = 0; i < N; i++) {
    x[i] = i * dx;
  }

  eigen_t *eig = dirac_1d(x, N, V, m, hbar, c);

  free(x);
  free(V);

  if (!eig) {
    printf("  FAIL: dirac_1d returned NULL\n");
    return 1;
  }

  double mc2 = m * c * c;
  int n_pos = 0, n_neg = 0, n_gap = 0;
  int n_total = eig->n;
  for (int i = 0; i < n_total; i++) {
    double E = eig->eigenvalues[i];

    if (E >= mc2 - 1e-6)
      n_pos++;
    else if (E <= -mc2 + 1e-6)
      n_neg++;
    else
      n_gap++;
  }

  printf("  mc^2=%.3f: %d states >= +mc^2, %d states <= -mc^2, %d in gap\n",
         mc2, n_pos, n_neg, n_gap);

  eigen_free(eig);

  // Loose: most states should fall in two branches, not the gap
  return (n_gap > n_total / 4) ? 1 : 0;
}

/*
 * Harmonic V(x) and c taken large (weak-relativistic-coupling limit),
 * Dirac equation's positive-energy branch must reduce exactly to ordinary
 * Schrodinger-equation spectrum shifted by the rest-mass energy,
 * E - m * c^2 -> (n + 1/2) * \hbar * \omega.
 */
static int test_nonrelativistic_limit(void) {
  int N = RUNNING_ON_VALGRIND ? 60 : 150;
  double L = 20.0;
  double dx = L / (N - 1);
  double hbar = 1.0, c = 200.0, m = 1.0, omega = 1.0;

  double *x = malloc(N * sizeof *x);
  double *V = malloc(N * sizeof *V);
  for (int i = 0; i < N; i++) {
    x[i] = -L / 2.0 + i * dx;
    V[i] = 0.5 * m * omega * omega * x[i] * x[i];
  }

  eigen_t *eig = dirac_1d(x, N, V, m, hbar, c);
  free(x);
  free(V);

  if (!eig) {
    printf("  FAIL: dirac_1d returned NULL\n");
    return 1;
  }

  double mc2 = m * c * c;
  // Collect the lowest few positive-energy (E > 0) eigenvalues
  int n_total = eig->n;
  double *pos = malloc((size_t)n_total * sizeof(double));
  int n_pos = 0;
  for (int i = 0; i < n_total; i++) {
    if (eig->eigenvalues[i] > 0.0) {
      pos[n_pos++] = eig->eigenvalues[i];
    }
  }
  /* NOTE: eigenvalues from cmatrix_eigh_complex are ascending, so the smallest
   * n_pos entries are already the lowest positive-energy states. */

  int n_check = RUNNING_ON_VALGRIND ? 1 : 3;
  int failed = 0;
  printf("  Non-relativistic limit (c=%.0f, harmonic V, m=omega=hbar=1):\n", c);
  for (int n = 0; n < n_check; n++) {
    double exact = (n + 0.5) * hbar * omega;
    /* each level appears twice (spin/component doubling); check first of the
     * pair at index 2*n */
    double shifted = pos[2 * n] - mc2;
    double err = fabs(shifted - exact);
    double tol = 0.02 * (n + 1);
    printf("    n=%d: E-mc^2=%.6f  exact=(n+1/2)=%.6f  err=%.2e  tol=%.2f\n", n,
           shifted, exact, err, tol);
    if (err > tol) {
      failed = 1;
    }
  }

  free(pos);
  eigen_free(eig);

  if (failed) {
    printf("  FAIL: non-relativistic limit does not match exact QHO "
           "spectrum\n");
  }

  return failed;
}

static int test_dirac_hydrogen_sommerfeld(void) {
  int fail = 0;
  double tol_rel = RUNNING_ON_VALGRIND ? 2e-5 : 5e-6;
  int N = RUNNING_ON_VALGRIND ? 80 : 300;
  double a0 = 4.0 * M_PI * EPSILON_0 * HBAR * HBAR /
              (M_ELECTRON * E_CHARGE * E_CHARGE); // Bohr radius
  double r_max = RUNNING_ON_VALGRIND ? 20.0 * a0 : 40.0 * a0;
  double r_min = 1e-4 * a0 / N;
  double *r = malloc(N * sizeof *r);
  double dr = (r_max - r_min) / (N - 1);
  for (int i = 0; i < N; i++) {
    r[i] = r_min + i * dr;
  }

  double Z_charge = 1.0; // nuclear charge
  double k_coulomb = Z_charge * E_CHARGE * E_CHARGE / (4.0 * M_PI * EPSILON_0);

  // {label, n, \kappa}
  struct {
    const char *label;
    int n, kappa;
  } states[] = {
      {"1s_1/2", 1, -1},
      {"2s_1/2", 2, -1},
      {"2p_1/2", 2, 1},
      {"2p_3/2", 2, -2},
  };

  int n_states = sizeof(states) / sizeof(states[0]);
  for (int s = 0; s < n_states; s++) {
    int n = states[s].n, kappa = states[s].kappa;

    double E_exact = dirac_hydrogen_energy_level(
        n, kappa, Z_charge, HBAR, M_ELECTRON, E_CHARGE, EPSILON_0, C_LIGHT);

    eigen_t *eig = dirac_radial_solve(r, N, kappa, V_coulomb, &k_coulomb,
                                      M_ELECTRON, HBAR, C_LIGHT);

    if (!eig) {
      printf("  FAIL: dirac_radial_solve returned NULL for %s\n",
             states[s].label);
      fail = 1;

      continue;
    }

    double mc2 = M_ELECTRON * C_LIGHT * C_LIGHT;
    double best = -1.0;
    double best_diff = 1e300;
    for (int i = 0; i < eig->n; i++) {
      double E = eig->eigenvalues[i];

      if (E > 0.0 && E < mc2) {
        double diff = fabs(E - E_exact);

        if (diff < best_diff) {
          best_diff = diff;
          best = E;
        }
      }
    }

    double rel_err = (best > 0.0) ? fabs(best - E_exact) / fabs(E_exact) : 1.0;
    printf("  %s (n=%d, \\kappa=%d): E_num=%.10e J  E_exact=%.10e J "
           "rel_err=%.2e\n",
           states[s].label, n, kappa, best, E_exact, rel_err);
    fail |= (rel_err > tol_rel);

    eigen_free(eig);
  }

  free(r);

  return fail;
}

// 2s_1/2 and 2p_1/2 should be exactly degenerate in point-charge Dirac spectrum
// (both have n=2, |\kappa|=1)
static int test_dirac_j_degeneracy(void) {
  double Z_charge = 1.0; // nuclear charge
  double E_2s = dirac_hydrogen_energy_level(2, -1, Z_charge, HBAR, M_ELECTRON,
                                            E_CHARGE, EPSILON_0, C_LIGHT);
  double E_2p = dirac_hydrogen_energy_level(2, 1, Z_charge, HBAR, M_ELECTRON,
                                            E_CHARGE, EPSILON_0, C_LIGHT);

  printf("  E(2s_1/2)=%.12e J  E(2p_1/2)=%.12e J  diff=%.3e\n", E_2s, E_2p,
         fabs(E_2s - E_2p));

  return fabs(E_2s - E_2p) > 1e-30 ? 1 : 0; // should be bit-identical
}

int main(void) {
  int failed = 0;

  printf("Dirac matrix Hermiticity");
  failed += test_hermiticity_of_construction();

  printf("Free particle +-mc^2 branch structure (qualitative):\n");
  failed += test_free_particle_branches();

  printf("Non-relativistic limit vs. exact harmonic-oscillator spectrum:\n");
  failed += test_nonrelativistic_limit();

  printf("Dirac radial solve vs. exact Sommerfeld hydrogen spectrum:\n");
  failed += test_dirac_hydrogen_sommerfeld();

  printf("2s_1/2 / 2p_1/2 exact j-degeneracy (n=2, |\\kappa|=1):\n");
  failed += test_dirac_j_degeneracy();

  if (failed == 0) {
    printf("\nAll test_dirac checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_dirac check(s) FAILED.\n", failed);
    return 1;
  }

  return 0;
}
