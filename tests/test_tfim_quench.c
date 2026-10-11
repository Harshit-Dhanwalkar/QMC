/*
 * Test: exact free-fermion quench of transverse-field Ising chain against exact
 * diagonalisation
 *
 * 1. t = 0 is ground state: <sx> equals ED ground state, L = 1, and
 *    correlations match ED
 * 2. After a sudden quench (two quenches, one across h = J in each direction)
 *    transverse magnetisation, Loschmidt echo and <sz_0 sz_r> for r = 1..3
 *    agree with an ED time evolution of same 8-site ring
 * 3. No quench: nothing moves
 * 4. Dynamical quantum phase transition: for a quench across h = J critical
 *    mode has overlap zero at t*_n, so rate function peaks there and critical
 *    time follows closed form
 * 5. Light cone: correlations beyond Lieb-Robinson front stay at their initial
 *    value
 * 6. Invalid input
 */

#include "../core/linalg/complex_eigh.h"
#include "../core/sparse.h"
#include "../physics/ising_chain.h"
#include "../physics/tfim_quench.h"
#include "complex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

enum { N = 8, DIM = 1 << N };

/* Even-parity ground state at field h on 8-site ring: lowest eigenvector from
 * ED, projected onto \sum_j sx_j parity +1 by flipping every spin (ordered
 * phase is nearly degenerate between two sectors) */
static int even_ground_state(double J, double h, complex_t *psi) {
  cmatrix_t *H = ising_hamiltonian_dense(N, J, h, 1);
  if (!H) {
    return -1;
  }

  eigen_t *eig = cmatrix_eigh_complex(H);
  if (!eig) {
    cmatrix_free(H);

    return -1;
  }

  double norm = 0.0;
  for (int s = 0; s < DIM; s++) {
    complex_t a = CMAT(eig->eigenvectors, s, 0),
              b = CMAT(eig->eigenvectors, (DIM - 1) ^ s, 0);

    psi[s] = c_add(a, b);
    norm += psi[s].re * psi[s].re + psi[s].im * psi[s].im;
  }

  for (int s = 0; s < DIM; s++) {
    psi[s].re /= sqrt(norm);
    psi[s].im /= sqrt(norm);
  }

  cmatrix_free(H);
  eigen_free(eig);

  return 0;
}

static double expect_sx(const complex_t *st) {
  double sum = 0.0;

  for (int j = 0; j < N; j++) {
    for (int s = 0; s < DIM; s++) {
      complex_t a = st[s], b = st[s ^ (1 << j)];

      sum += a.re * b.re + a.im * b.im;
    }
  }

  return sum / N;
}

static double expect_zz(const complex_t *st, int r) {
  double sum = 0.0;

  for (int s = 0; s < DIM; s++) {
    double z0 = (s & 1) ? -1.0 : 1.0, zr = ((s >> r) & 1) ? -1.0 : 1.0;

    sum += z0 * zr * (st[s].re * st[s].re + st[s].im * st[s].im);
  }

  return sum;
}

static void compare_with_ed(double J, double h_i, double h_f) {
  printf("  --- quench h = %.2f -> %.2f against ED (%d sites) ---\n", h_i, h_f,
         N);

  static complex_t psi[DIM], st[DIM], coef[DIM];

  if (even_ground_state(J, h_i, psi) != 0) {
    check_true(0, "ground state");

    return;
  }

  cmatrix_t *Hf = ising_hamiltonian_dense(N, J, h_f, 1);
  eigen_t *eig = Hf ? cmatrix_eigh_complex(Hf) : NULL;
  if (!eig) {
    check_true(0, "post-quench diagonalisation");
    cmatrix_free(Hf);

    return;
  }

  for (int e = 0; e < DIM; e++) {
    complex_t sum = {0.0, 0.0};

    for (int s = 0; s < DIM; s++) {
      sum = c_add(sum, c_mul(c_conj(CMAT(eig->eigenvectors, s, e)), psi[s]));
    }

    coef[e] = sum;
  }

  const double times[] = {0.0, 0.3, 0.9, 1.6};

  for (int it = 0; it < 4; it++) {
    double t = times[it];

    for (int s = 0; s < DIM; s++) {
      complex_t sum = {0.0, 0.0};

      for (int e = 0; e < DIM; e++) {
        double ph = -eig->eigenvalues[e] * t;
        complex_t u = {cos(ph), sin(ph)};

        sum =
            c_add(sum, c_mul(CMAT(eig->eigenvectors, s, e), c_mul(u, coef[e])));
      }

      st[s] = sum;
    }

    complex_t ov = {0.0, 0.0};

    for (int s = 0; s < DIM; s++) {
      ov = c_add(ov, c_mul(c_conj(psi[s]), st[s]));
    }

    double L = ov.re * ov.re + ov.im * ov.im;
    double zz[4];
    char label[96];

    tfim_quench_zz(J, h_i, h_f, t, N, 3, zz);
    snprintf(label, sizeof label, "t=%.1f <sx>", t);
    check_close(tfim_quench_mx(J, h_i, h_f, t, N), expect_sx(st), 1e-9, label);

    snprintf(label, sizeof label, "t=%.1f rate function", t);
    check_close(tfim_quench_loschmidt_rate(J, h_i, h_f, t, N), -log(L) / N,
                1e-8, label);

    for (int r = 1; r <= 3; r++) {
      snprintf(label, sizeof label, "t=%.1f <sz0 sz%d>", t, r);
      check_close(zz[r], expect_zz(st, r), 1e-9, label);
    }
  }

  cmatrix_free(Hf);
  eigen_free(eig);
}

int main(void) {
  printf(" > TTFIM Quench tests\n");

  const double J = 1.0;

  compare_with_ed(J, 0.3, 1.7);
  compare_with_ed(J, 1.8, 0.4);

  printf("  --- no quench: nothing moves ---");
  double zz0[8];
  double zz1[8];

  tfim_quench_zz(J, 0.8, 0.8, 0.0, 64, 4, zz0);
  tfim_quench_zz(J, 0.8, 0.8, 7.3, 64, 4, zz1);
  check_close(tfim_quench_mx(J, 0.8, 0.8, 5.0, 64),
              tfim_quench_mx(J, 0.8, 0.8, 0.0, 64), 1e-12, "<sx> constant");
  check_close(zz1[3], zz0[3], 1e-12, "<sz0 sz3> constant");
  check_close(tfim_quench_loschmidt_rate(J, 0.8, 0.8, 3.0, 64), 0.0, 1e-12,
              "rate function zero");

  printf("  --- deep paramagnet / ferromagnet limits ---");
  check_close(tfim_quench_mx(J, 1e9, 1e9, 0.0, 32), 1.0, 1e-6,
              "<sx> -> 1 for h -> infinity");
  tfim_quench_zz(J, 1e-9, 1e-9, 0.0, 32, 5, zz0);
  check_close(zz0[5], 1.0, 1e-6, "<sz0 sz5> -> 1 for h -> 0");

  printf("  --- dynamical quantum phase transition (h = 0.3 -> 2.0) ---");
  {
    double hi = 0.3, hf = 2.0, t0 = tfim_quench_critical_time(J, hi, hf, 0);
    double t1 = tfim_quench_critical_time(J, hi, hf, 1);
    double c = (J * J + hi * hf) / (J * (hi + hf)), kstar = acos(c);
    double eps = 2.0 * sqrt(J * J + hf * hf - 2.0 * J * hf * cos(kstar));
    int n_big = 4096;

    check_close(t0, M_PI / (2.0 * eps), 1e-12, "t*_0 = \\pi / 2 \\eps_k*");
    check_close(t1 / t0, 3.0, 1e-12, "t*_1 / t*_0 = 3");

    /* rate function is smooth on either side but has a kink at t*: its slope
     * changes sign across it, so check maximum of a fine scan sits at t*_0
     * within scan step */
    double best_t = 0.0, best = -1.0, step = 0.002;

    for (double t = 0.05; t < 1.5 * t0; t += step) {
      double lam = tfim_quench_loschmidt_rate(J, hi, hf, t, n_big);

      if (lam > best) {
        best = lam;
        best_t = t;
      }
    }

    printf("  peak of rate function at t=%.4f, t*_0=%.4f\n", best_t, t0);
    check_true(fabs(best_t - t0) < 0.02,
               "rate function peaks at critical time");
    check_true(isnan(tfim_quench_critical_time(J, 0.3, 0.6, 0)),
               "no critical time without crossing h = J");
  }

  printf("  --- light cone ---");
  {
    /* start deep in paramagnet, where <sz_0 sz_r> ~ (J/h_i)^r is
     * negligible, and quench into ordered phase: correlations spread
     * ballistically at twice the maximum quasiparticle speed 2 min(J, h_f) */
    double hi = 3.0, hf = 0.5, t = 3.0, zz[41];
    int n = 256;

    tfim_quench_zz(J, hi, hf, t, n, 40, zz);
    double speed = 4.0 * (hf < J ? hf : J), front = speed * t;

    printf("  front at r = %.1f\n", front);
    check_true(fabs(zz[3]) > 0.05, "correlations have changed inside the cone");
    for (int r = (int)(front + 6.0); r <= 40; r += 8) {
      char label[64];

      snprintf(label, sizeof label, "r=%d outside the cone is ~0", r);
      check_close(zz[r], 0.0, 2e-3, label);
    }
  }

  printf("  --- invalid input ---");
  double dummy[10];

  check_true(isnan(tfim_quench_mx(J, 0.5, 1.0, 0.0, 7)),
             "odd chain length rejected");
  check_true(isnan(tfim_quench_mx(J, NAN, 1.0, 0.0, 8)), "NaN field rejected");
  check_true(tfim_quench_zz(J, 0.5, 1.0, 0.0, 8, 4, dummy) == -1,
             "r_max beyond half the ring rejected");
  check_true(tfim_quench_zz(J, 0.5, 1.0, 0.0, 64, 70, dummy) == -1,
             "r_max beyond the limit rejected");
  check_true(tfim_quench_zz(J, 0.5, 1.0, 0.0, 64, 3, NULL) == -1,
             "NULL output rejected");
  check_true(tfim_quench_modes(J, 0.5, 1.0, 0.0, 8, NULL) == -1,
             "NULL modes rejected");

  if (failures) {
    printf("\n  %d check(s) FAILED\n", failures);

    return 1;
  }
  printf("\n  all checks passed\n");

  return 0;
}
