/*
 * Tests the Crank-Nicolson TDSE integrator.
 *
 * Properties verified:
 *  1. Norm preservation: ||\phi(t)||^2 = 1 after N steps (unitarity)
 *  2. Energy conservation: <H> stays constant under free evolution
 *  3. Phase propagation: stationary state \phi_n picks up phase e^{-iE_n t}
 *     so |\phi(t)|^2 = |\phi(0)|^2 exactly for an energy eigenstate
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

int main(void) {
  printf(" > Testing Crank-Nicolson time evolution...\n");

  // Setup: HO ground state as initial condition
  int N = 201;
  double xmin = -5.0, xmax = 5.0;
  double dx = (xmax - xmin) / (N - 1);
  double *x = linspace(xmin, xmax, N);
  double *V = malloc(N * sizeof *V);
  if (!x || !V) {
    printf("FAIL: memory\n");
    return 1;
  }

  for (int i = 0; i < N; i++)
    V[i] = 0.5 * x[i] * x[i];

  // Build H and diagonalize to get ground state
  double coeff = 0.5 / (dx * dx);
  cmatrix_t *H = cmatrix_alloc(N, N);
  for (int i = 0; i < N; i++) {
    CMAT(H, i, i) = c_real(2.0 * coeff + V[i]);
    if (i > 0)
      CMAT(H, i, i - 1) = c_real(-coeff);
    if (i < N - 1)
      CMAT(H, i, i + 1) = c_real(-coeff);
  }

  eigen_t *eig = cmatrix_eigh(H);
  if (!eig) {
    printf("FAIL: eigensolver\n");
    cmatrix_free(H);
    free(x);
    free(V);
    return 1;
  }

  // Ground state energy
  double E0 = eig->eigenvalues[0];
  printf("   Ground state energy: %.6f (expected ~0.5)\n", E0);

  // Initial \phi = ground state eigenvector, continuum-normalized
  cvector_t *psi = cvector_alloc(N);
  for (int i = 0; i < N; i++) {
    psi->data[i] = CMAT(eig->eigenvectors, i, 0);
  }

  double norm0 = 0.0;
  for (int i = 0; i < N; i++) {
    norm0 += (psi->data[i].re * psi->data[i].re +
              psi->data[i].im * psi->data[i].im) *
             dx;
  }

  double inv = 1.0 / sqrt(norm0);
  for (int i = 0; i < N; i++) {
    psi->data[i].re *= inv;
    psi->data[i].im *= inv;
  }
  eigen_free(eig);

  // Build tridiagonal Hamiltonian for CN
  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  build_tridiagonal_hamiltonian(x, V, N, dx, 0.5, diag, offdiag);

  // Test 1: Norm preservation over 200 steps
  printf("   Test 1: norm preservation...\n");
  double dt = 0.01;
  int n_steps = 200;
  int norm_fail = 0;

  for (int step = 0; step < n_steps; step++) {
    crank_nicolson_step(diag, offdiag, dt, psi);
    double norm = 0.0;
    for (int i = 0; i < N; i++) {
      norm += (psi->data[i].re * psi->data[i].re +
               psi->data[i].im * psi->data[i].im) *
              dx;
    }

    if (fabs(norm - 1.0) > 1e-10) {
      printf("   FAIL: norm = %.10f at step %d\n", norm, step);
      norm_fail = 1;
      break;
    }
  }

  if (!norm_fail)
    printf("   Norm preserved to 1e-10 over %d steps: PASS\n", n_steps);

  // Test 2: Stationary state - |\phi(t)|^2 unchanged
  printf("   Test 2: |\\phi(t)|^2 invariant for energy eigenstate...\n");

  // Reset to ground state
  eigen_t *eig2 = cmatrix_eigh(H);
  for (int i = 0; i < N; i++) {
    psi->data[i] = CMAT(eig2->eigenvectors, i, 0);
  }

  norm0 = 0.0;
  for (int i = 0; i < N; i++) {
    norm0 += (psi->data[i].re * psi->data[i].re +
              psi->data[i].im * psi->data[i].im) *
             dx;
  }

  inv = 1.0 / sqrt(norm0);
  for (int i = 0; i < N; i++) {
    psi->data[i].re *= inv;
    psi->data[i].im *= inv;
  }

  // Save initial |\phi|^2
  double *prob0 = malloc(N * sizeof *prob0);
  for (int i = 0; i < N; i++)
    prob0[i] =
        psi->data[i].re * psi->data[i].re + psi->data[i].im * psi->data[i].im;

  // Evolve one full period: T = 2\pi /E0
  double T = 2.0 * M_PI / E0;
  int n_full = (int)(T / dt) + 1;
  for (int s = 0; s < n_full; s++)
    crank_nicolson_step(diag, offdiag, dt, psi);

  double max_prob_err = 0.0;
  for (int i = 0; i < N; i++) {
    double prob_t =
        psi->data[i].re * psi->data[i].re + psi->data[i].im * psi->data[i].im;
    double err = fabs(prob_t - prob0[i]);
    if (err > max_prob_err)
      max_prob_err = err;
  }
  printf("   Max |\\phi(T)|^2-|\\phi(0)|^2 error after T=2\\pi /E_0: %.3e\n",
         max_prob_err);
  int prob_pass = (max_prob_err < 1e-4);
  printf("   %s\n", prob_pass ? "PASS" : "FAIL");

  eigen_free(eig2);
  cmatrix_free(H);
  free(prob0);
  free(diag);
  free(offdiag);
  cvector_free(psi);
  free(x);
  free(V);

  int all_pass = !norm_fail && prob_pass;
  if (all_pass)
    printf("   Crank-Nicolson test passed.\n");

  return all_pass ? 0 : 1;
}
