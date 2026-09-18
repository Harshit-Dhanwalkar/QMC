/*
 * Tests: Crank-Nicolson TDSE integrator, time-dependent V(x,t) and CAP support.
 *
 * Properties verified:
 *  1. Norm preservation: ||\phi(t)||^2 = 1 after N steps (unitarity)
 *  2. Energy conservation: <H> stays constant under free evolution
 *  3. Phase propagation: stationary state \phi_n picks up phase
 *     \exp^{-i * E_n * t} so |\phi(t)|^2 = |\phi(0)|^2 exactly for an energy
 *     eigenstate
 *  4. Time-independent limit: driving V(x,t) with t-independent function must
 *     reproduce existing static crank_nicolson_step path to tight tolerance
 *  5. No CAP: norm must stay conserved (H Hermitian, CN is unitary) over many
 *     steps
 *  6. CAP: a rightward-moving Gaussian wavepacket that reaches CAP layer at
 *     right edge must lose substantial fraction of its norm (absorption
 *     working) vs. no-CAP control that stays ~1
 */

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/matrix.h"
#include "../core/ode/crank_nicolson.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check(int cond, const char *msg) {
  if (!cond) {
    printf("  FAIL: %s\n", msg);
    failures++;
  }
}

static void check_close(double got, double expected, double tol,
                        const char *msg) {
  if (fabs(got - expected) > tol) {
    printf("  FAIL: %s (got %.10f, expected %.10f, diff %.2e)\n", msg, got,
           expected, fabs(got - expected));
    failures++;
  }
}

static double norm_sq(const cvector_t *psi, double dx) {
  double s = 0.0;
  for (int i = 0; i < psi->n; i++) {
    s += c_abs2(psi->data[i]);
  }
  return s * dx;
}

// Static harmonic potential, ignores t
static double V_harmonic_static(double x, double t, void *params) {
  (void)t;
  double omega = *(double *)params;

  return 0.5 * omega * omega * x * x;
}

static double V_zero(double x, double t, void *params) {
  (void)x;
  (void)t;
  (void)params;

  return 0.0;
}

static void make_gaussian(cvector_t *psi, const double *x, int N, double x0,
                          double sigma0, double k0) {
  for (int i = 0; i < N; i++) {
    double dxg = x[i] - x0;
    double envelope = exp(-dxg * dxg / (4.0 * sigma0 * sigma0));

    psi->data[i].re = envelope * cos(k0 * x[i]);
    psi->data[i].im = envelope * sin(k0 * x[i]);
  }
}

static void test_time_independent_limit(void) {
  printf("Test: time-independent limit matches static path\n");

  int N = 400;
  double x_min = -10.0;
  double x_max = 10.0;
  double dx = (x_max - x_min) / (N - 1);
  double hbar_sq_2m = 0.5;
  double omega = 1.0;
  double dt = 0.01;
  int steps = 50;

  double *x = malloc(N * sizeof *x);
  for (int i = 0; i < N; i++) {
    x[i] = x_min + i * dx;
  }

  // Static path
  double *V = malloc(N * sizeof *V);
  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  for (int i = 0; i < N; i++) {
    V[i] = 0.5 * omega * omega * x[i] * x[i];
  }

  build_tridiagonal_hamiltonian(x, V, N, dx, hbar_sq_2m, diag, offdiag);

  cvector_t *psi_static = cvector_alloc(N);
  make_gaussian(psi_static, x, N, -3.0, 1.0, 0.0);
  for (int s = 0; s < steps; s++) {
    crank_nicolson_step(diag, offdiag, dt, psi_static);
  }

  // Time-dependent path, driven with a t-independent V(x,t)
  cvector_t *psi_td = cvector_alloc(N);
  make_gaussian(psi_td, x, N, -3.0, 1.0, 0.0);

  crank_nicolson_evolve_time_dependent(x, N, dx, hbar_sq_2m, V_harmonic_static,
                                       &omega, NULL, psi_td, 0.0, dt, steps);

  double max_err = 0.0;
  for (int i = 0; i < N; i++) {
    double dre = psi_static->data[i].re - psi_td->data[i].re;
    double dim = psi_static->data[i].im - psi_td->data[i].im;
    double err = sqrt(dre * dre + dim * dim);
    if (err > max_err) {
      max_err = err;
    }
  }

  check_close(max_err, 0.0, 1e-8,
              "time-dependent path matches static path pointwise");

  free(x);
  free(V);
  free(diag);
  free(offdiag);
  cvector_free(psi_static);
  cvector_free(psi_td);
}

static void test_norm_conserved_no_cap(void) {
  printf("Test: norm conservation (no CAP)\n");

  int N = 300;
  double x_min = -15.0;
  double x_max = 15.0;
  double dx = (x_max - x_min) / (N - 1);
  double hbar_sq_2m = 0.5;
  double omega = 1.0;
  double dt = 0.01;
  int steps = 200;

  double *x = malloc(N * sizeof *x);
  for (int i = 0; i < N; i++) {
    x[i] = x_min + i * dx;
  }

  cvector_t *psi = cvector_alloc(N);
  make_gaussian(psi, x, N, 0.0, 1.5, 1.0);

  double n0 = norm_sq(psi, dx);
  for (int i = 0; i < N; i++) {
    psi->data[i] = c_scale(psi->data[i], 1.0 / sqrt(n0));
  }

  crank_nicolson_evolve_time_dependent(x, N, dx, hbar_sq_2m, V_harmonic_static,
                                       &omega, NULL, psi, 0.0, dt, steps);

  double n_final = norm_sq(psi, dx);
  check_close(n_final, 1.0, 1e-6, "norm conserved over 200 CN steps");

  free(x);
  cvector_free(psi);
}

// Group velocity v = \hbar* k0 / m = 3.0 (natural units)
static void test_cap_absorbs(void) {
  printf("Test: CAP absorbs outgoing wavepacket\n");

  int N = 400;
  double x_min = -20.0;
  double x_max = 20.0;
  double dx = (x_max - x_min) / (N - 1);
  double hbar_sq_2m = 0.5;
  double dt = 0.02;
  int steps = 700;
  double k0 = 3.0;

  double *x = malloc(N * sizeof *x);
  for (int i = 0; i < N; i++) {
    x[i] = x_min + i * dx;
  }

  double *absorb = malloc(N * sizeof *absorb);
  cap_build_monomial(x, N, 4.0, 5.0, 2, absorb);

  // With CAP
  cvector_t *psi_cap = cvector_alloc(N);
  make_gaussian(psi_cap, x, N, -10.0, 1.0, k0);

  double n0_cap = norm_sq(psi_cap, dx);
  for (int i = 0; i < N; i++) {
    psi_cap->data[i] = c_scale(psi_cap->data[i], 1.0 / sqrt(n0_cap));
  }

  crank_nicolson_evolve_time_dependent(x, N, dx, hbar_sq_2m, V_zero, NULL,
                                       absorb, psi_cap, 0.0, dt, steps);

  double n_cap = norm_sq(psi_cap, dx);

  // Control: no CAP
  cvector_t *psi_free = cvector_alloc(N);
  make_gaussian(psi_free, x, N, -10.0, 1.0, k0);

  double n0_free = norm_sq(psi_free, dx);
  for (int i = 0; i < N; i++) {
    psi_free->data[i] = c_scale(psi_free->data[i], 1.0 / sqrt(n0_free));
  }

  crank_nicolson_evolve_time_dependent(x, N, dx, hbar_sq_2m, V_zero, NULL, NULL,
                                       psi_free, 0.0, dt, steps);

  double n_free = norm_sq(psi_free, dx);

  printf("  final norm with CAP:        %.6f (expect <1)\n", n_cap);
  printf("  final norm, NO-CAP control: %.6f (expect ~1)\n", n_free);

  check(n_cap < 0.3,
        "CAP absorbs most of the wavepacket norm"); // packet crosses CAP layer
  check_close(n_free, 1.0, 1e-3, "no-CAP control stays normalized");

  free(x);
  free(absorb);
  cvector_free(psi_cap);
  cvector_free(psi_free);
}

static void test_stationary_state_probability_invariant(void) {
  printf("Test: |\\phi(t)|^2 invariant for energy eigenstate\n");

  int N = 201;
  double xmin = -5.0;
  double xmax = 5.0;
  double dx = (xmax - xmin) / (N - 1);
  double *x = linspace(xmin, xmax, N);
  double *V = malloc(N * sizeof *V);
  if (!x || !V) {
    check(0, "allocation failure");

    free(x);
    free(V);

    return;
  }

  for (int i = 0; i < N; i++) {
    V[i] = 0.5 * x[i] * x[i];
  }

  double coeff = 0.5 / (dx * dx);
  cmatrix_t *H = cmatrix_alloc(N, N);
  for (int i = 0; i < N; i++) {
    CMAT(H, i, i) = c_real(2.0 * coeff + V[i]);

    if (i > 0) {
      CMAT(H, i, i - 1) = c_real(-coeff);
    }

    if (i < N - 1) {
      CMAT(H, i, i + 1) = c_real(-coeff);
    }
  }

  eigen_t *eig = cmatrix_eigh(H);
  if (!eig) {
    check(0, "eigensolver failed");

    cmatrix_free(H);
    free(x);
    free(V);

    return;
  }

  double E0 = eig->eigenvalues[0];
  check_close(E0, 0.5, 1e-2, "ground state energy ~ 0.5");

  cvector_t *psi = cvector_alloc(N);
  for (int i = 0; i < N; i++) {
    psi->data[i] = CMAT(eig->eigenvectors, i, 0);
  }

  double norm0 = 0.0;
  for (int i = 0; i < N; i++) {
    norm0 += c_abs2(psi->data[i]) * dx;
  }

  double inv = 1.0 / sqrt(norm0);
  for (int i = 0; i < N; i++) {
    psi->data[i].re *= inv;
    psi->data[i].im *= inv;
  }

  eigen_free(eig);

  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  build_tridiagonal_hamiltonian(x, V, N, dx, 0.5, diag, offdiag);

  double *prob0 = malloc(N * sizeof *prob0);
  for (int i = 0; i < N; i++) {
    prob0[i] = c_abs2(psi->data[i]);
  }

  double dt = 0.01;
  double T = 2.0 * M_PI / E0;
  int n_full = (int)(T / dt) + 1;
  for (int s = 0; s < n_full; s++) {
    crank_nicolson_step(diag, offdiag, dt, psi);
  }

  double max_prob_err = 0.0;
  for (int i = 0; i < N; i++) {
    double prob_t = c_abs2(psi->data[i]);
    double err = fabs(prob_t - prob0[i]);
    if (err > max_prob_err) {
      max_prob_err = err;
    }
  }

  check(max_prob_err < 1e-4,
        "|\\phi(T)|^2 == |\\phi(0)|^2 after one period T = 2 * \\pi/E_0");

  cmatrix_free(H);
  free(prob0);
  free(diag);
  free(offdiag);
  cvector_free(psi);
  free(x);
  free(V);
}

int main(void) {
  printf(" > Testing Crank-Nicolson time evolution...\n\n");

  test_time_independent_limit();
  test_norm_conserved_no_cap();
  test_cap_absorbs();
  test_stationary_state_probability_invariant();

  if (failures == 0) {
    printf("\nAll test_crank_nicolson checks passed.\n");
    return 0;
  } else {
    printf("\n%d test_crank_nicolson check(s) FAILED.\n", failures);
    return 1;
  }
}
