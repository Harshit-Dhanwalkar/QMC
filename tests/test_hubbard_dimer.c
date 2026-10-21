/*
 * Test: Hubbard dimer by exact diagonalisation
 *
 * 1. N = 2 spectrum: U/2 -+ \sqrt(U^2/4 + 4 t^2), U and a 3-fold zero;
 *    trace and eigenvector orthonormality
 * 2. Ground state: ionic weight (1 - U/2\Delta)/2, <S1.S2> = -(3/4)(1 - ionic),
 *    one-site entropy from 4 probabilities
 * 3. Limits: U = 0 gives <S1.S2> = -3/8, large U gives Heisenberg pair (E0 ->
 *    -4 t^2/U, <S1.S2> -> -3/4)
 * 4. Eigenstates have definite N and S_z; triplet has <S1.S2> = 1/4 and
 *    <S^2>-type consistency
 * 5. Quench from Neel state: norm and energy conserved, staggered moment
 *    oscillates as \cos(J t) for U >> t and long-time doublon number matches
 *    diagonal ensemble
 * 6. Canonical thermal averages: T -> 0 is ground state, T -> \infty  gives D =
 *    1/3, <S1.S2> = 0, S = \ln(6)
 * 7. Invalid input
 */

#include "../physics/hubbard_dimer.h"
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int failures = 0;

static void check_close(double got, double expected, double tol,
                        const char *label) {
  double err = fabs(got - expected);

  printf("  %s: got=%.8g expected=%.8g err=%.2e\n", label, got, expected, err);
  if (!(err <= tol)) {
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

static double norm2(const double *re, const double *im) {
  double s = 0.0;

  for (int m = 0; m < HD_DIM; m++) {
    s += re[m] * re[m] + im[m] * im[m];
  }

  return s;
}

// Main test driver for Hubbard dimer exact-diagonalisation routines
//
// Checks below validate:
// - analytical N=2 spectrum at representative U values
// - global Hamiltonian invariants (trace in full Fock space)
// - eigenvector orthonormality from numerical solver
// - additional observable/state tests in later sections of this function
int main(void) {
  const double us[] = {0.0, 1.0, 4.0, 10.0};

  // Compare computed N=2 energies against known closed-form spectrum
  printf("  === Test: Spectrum ===\n");
  for (int c = 0; c < 4; c++) {
    double u = us[c];
    double t = 1.0;
    double d = sqrt(u * u / 4 + 4 * t * t);
    hd_solution_t s;
    int idx[6];

    hd_solve(t, u, &s);
    hd_n2_states(&s, idx);
    double want[6] = {u / 2 - d, 0, 0, 0, u, u / 2 + d};

    // Sort they can be compared index-by-index
    for (int i = 0; i < 6; i++) {
      for (int j = i + 1; j < 6; j++) {
        if (want[j] < want[i]) {
          double tmp = want[i];

          want[i] = want[j];
          want[j] = tmp;
        }
      }
    }

    for (int i = 0; i < 6; i++) {
      char lab[64];

      snprintf(lab, sizeof lab, "U=%g level %d", u, i);
      check_close(s.e[idx[i]], want[i], 1e-10, lab);
    }

    // Trace(H) equals U times number of doublon-containing basis states
    double tr = 0.0;

    for (int k = 0; k < HD_DIM; k++) {
      tr += s.e[k];
    }

    check_close(tr, 8 * u, 1e-9,
                "trace of H over Fock space (8 configs hold a doublon)");
  }

  // Verify that eigenvectors returned by solver form an orthonormal basis
  {
    hd_solution_t s;
    double worst = 0.0;

    hd_solve(1.0, 3.0, &s);
    for (int a = 0; a < HD_DIM; a++) {
      for (int b = 0; b < HD_DIM; b++) {
        double dot = 0.0;

        for (int m = 0; m < HD_DIM; m++) {
          dot += s.v[a][m] * s.v[b][m];
        }

        worst = fmax(worst, fabs(dot - (a == b)));
      }
    }

    check_true(worst < 1e-12, "eigenvectors orthonormal");
  }

  printf("  === Test: Ground state closed forms ===\n");
  for (int c = 0; c < 4; c++) {
    double u = us[c];
    double t = 1.0;
    double d = sqrt(u * u / 4 + 4 * t * t);
    hd_solution_t s;
    int idx[6];
    double re[HD_DIM];
    double im[HD_DIM];
    double p[4];

    hd_solve(t, u, &s);
    hd_n2_states(&s, idx);
    hd_eigenstate(&s, idx[0], re, im);

    double ion = 0.5 * (1.0 - u / (2 * d));

    check_close(hd_expect_double_occ(re, im), ion, 1e-10, "double occupancy");
    check_close(hd_expect_ss(re, im), -0.75 * (1.0 - ion), 1e-10, "<S1.S2>");
    check_close(hd_expect_sz(re, im, 0), 0.0, 1e-12, "<Sz1>");

    hd_site_probs(re, im, 0, p);

    check_close(p[0], ion / 2, 1e-10, "p(empty)");
    check_close(p[3], ion / 2, 1e-10, "p(double)");
    check_close(p[1], (1 - ion) / 2, 1e-10, "p(up)");

    double S = 0.0;

    for (int i = 0; i < 4; i++) {
      if (p[i] > 0) {
        S -= p[i] * log(p[i]);
      }
    }

    check_close(hd_site_entropy(re, im, 0), S, 1e-12, "site entropy");
    check_close(hd_site_entropy(re, im, 1), S, 1e-10, "both sites equal");
    check_close(hd_expect_energy(&s, re, im), u / 2 - d, 1e-10, "<H>");
  }

  printf("  === Test: Limits ===\n");
  {
    hd_solution_t s;
    int idx[6];
    double re[HD_DIM];
    double im[HD_DIM];

    hd_solve(1.0, 0.0, &s);
    hd_n2_states(&s, idx);
    hd_eigenstate(&s, idx[0], re, im);

    check_close(hd_expect_ss(re, im), -0.375, 1e-10, "U=0 <S1.S2>");
    check_close(hd_expect_double_occ(re, im), 0.5, 1e-10, "U=0 doublons");

    hd_solve(1.0, 200.0, &s);
    hd_n2_states(&s, idx);
    hd_eigenstate(&s, idx[0], re, im);

    check_close(s.e[idx[0]] * 200.0 / 4.0, -1.0, 1e-3, "E0 U / 4t^2 -> -1");
    check_close(hd_expect_ss(re, im), -0.75, 1e-3, "large-U <S1.S2>");
    check_close(hd_site_entropy(re, im, 0), log(2.0), 5e-3,
                "large-U entropy ln 2");
  }

  printf("  === Test: Definite N and Sz, triplet ===\n");
  {
    hd_solution_t s;
    int idx[6];
    double trip = 0;

    hd_solve(1.0, 4.0, &s);
    hd_n2_states(&s, idx);
    for (int i = 0; i < 6; i++) {
      double re[HD_DIM];
      double im[HD_DIM];

      hd_eigenstate(&s, idx[i], re, im);
      if (fabs(s.e[idx[i]]) < 1e-9) {
        trip++;

        check_close(hd_expect_ss(re, im), 0.25, 1e-10, "triplet <S1.S2>");
      }
    }

    check_true(trip == 3, "three triplet states");
    check_true(s.n[0] == 1 && fabs(s.e[0] + 1.0) < 1e-10,
               "lowest Fock-space level at U = 4 is a single electron at -t");
  }

  printf("  === Test: Neel quench ===\n");
  {
    hd_solution_t s;
    double re0[HD_DIM];
    double im0[HD_DIM];
    double re[HD_DIM];
    double im[HD_DIM];
    double e0;
    const double u = 40.0;
    double j = 4.0 / u;

    hd_solve(1.0, u, &s);
    hd_basis_state(9, re0, im0); // 1 up, 2 down

    e0 = hd_expect_energy(&s, re0, im0);

    check_close(hd_expect_sz(re0, im0, 0), 0.5, 1e-12, "Neel Sz1");
    check_close(hd_expect_sz(re0, im0, 1), -0.5, 1e-12, "Neel Sz2");

    double worst_m = 0.0;

    for (int i = 1; i <= 12; i++) {
      double time = 0.5 * M_PI / j * i / 12.0 * 4.0;

      hd_evolve(&s, re0, im0, time, re, im);
      if (fabs(norm2(re, im) - 1.0) > 1e-10 ||
          fabs(hd_expect_energy(&s, re, im) - e0) > 1e-9) {
        failures++;
        printf("  FAIL: norm/energy drift at t=%g\n", time);
      }

      worst_m =
          fmax(worst_m, fabs(2 * hd_expect_sz(re, im, 0) - cos(j * time)));
    }

    printf("  max |2<Sz1>(t) - \\cos(Jt)| over 2 periods-ish = %.2e\n",
           worst_m);
    check_true(worst_m < 0.02, "staggered moment follows \\cos(Jt)");
    hd_evolve(&s, re0, im0, 0.0, re, im);

    check_close(hd_expect_double_occ(re, im), 0.0, 1e-12,
                "Neel has no doublons at t=0");

    // diagonal ensemble of the doublon number
    double diag = 0.0;
    int idx[6];

    hd_n2_states(&s, idx);
    for (int i = 0; i < 6; i++) {
      double rk[HD_DIM];
      double ik[HD_DIM];
      double ov = 0.0;

      hd_eigenstate(&s, idx[i], rk, ik);
      for (int m = 0; m < HD_DIM; m++) {
        ov += rk[m] * re0[m];
      }

      diag += ov * ov * hd_expect_double_occ(rk, ik);
    }

    double avg = 0.0;
    const int steps = 40000;

    for (int i = 0; i < steps; i++) {
      hd_evolve(&s, re0, im0, 0.37 + i * 0.731, re, im);
      avg += hd_expect_double_occ(re, im) / steps;
    }

    check_close(avg, diag, 5e-4, "time-averaged doublons = diagonal ensemble");
  }

  printf("  === Test: Thermal averages ===\n");
  {
    hd_solution_t s;
    hd_thermal_t th;
    int idx[6];
    double re[HD_DIM];
    double im[HD_DIM];

    hd_solve(1.0, 4.0, &s);
    hd_n2_states(&s, idx);
    hd_eigenstate(&s, idx[0], re, im);
    hd_thermal_n2(&s, 0.01, &th);

    check_close(th.ss, hd_expect_ss(re, im), 1e-9, "T->0 <S1.S2>");
    check_close(th.double_occ, hd_expect_double_occ(re, im), 1e-9,
                "T->0 doublons");
    check_close(th.entropy, 0.0, 1e-9, "T->0 entropy");

    hd_thermal_n2(&s, 1e6, &th);

    check_close(th.double_occ, 1.0 / 3.0, 1e-4, "T->inf doublons 1/3");
    check_close(th.ss, 0.0, 1e-4, "T->inf <S1.S2>");
    check_close(th.entropy, log(6.0), 1e-4, "T->inf entropy ln 6");

    hd_thermal_n2(&s, 0.3, &th);

    hd_thermal_t hot;

    hd_thermal_n2(&s, 0.6, &hot);

    check_true(hot.energy > th.energy && hot.entropy > th.entropy &&
                   hot.ss > th.ss,
               "energy, entropy, <S1.S2> rise with T");
  }

  printf("  === Test: Invalid input ===\n");
  {
    hd_solution_t s;
    hd_thermal_t th;
    double re[HD_DIM];
    double im[HD_DIM];

    check_true(hd_solve(NAN, 1.0, &s) == -1 &&
                   hd_solve(1.0, INFINITY, &s) == -1,
               "non-finite parameters rejected");

    hd_solve(1.0, 2.0, &s);
    hd_basis_state(9, re, im);

    check_true(hd_evolve(&s, re, im, NAN, re, im) == -1, "NaN time rejected");
    check_true(hd_thermal_n2(&s, 0.0, &th) == -1 &&
                   hd_thermal_n2(&s, -1.0, &th) == -1,
               "non-positive temperature rejected");
  }

  printf(failures ? "\n  %d FAILED\n" : "\n  all checks passed\n", failures);

  return failures != 0;
}
