/*
Test + demonstration: Lindblad equation (open quantum systems).

1. Amplitude damping (T1) of a single excited qubit, H=0: population
   \rho11(t) decays as \exp(-\gamma * t) with no unitary drive - exact
   closed-form solution of d(\rho11)/dt = -\gamma * \rho11, lindblad_evolve's
   RK4 integration.
2. Pure dephasing (T2) of a |+> superposition, H=0: coherence \rho01(t)
   decays as \rho01(0) * exp(-\gamma * t) while populations stay fixed at 0.5
   each (derived analytically from the Lindblad dissipator for L =
   \sqrt(\gamma/2) * \sigma_z)
3. Closed-system cross-check (n_ops=0): density-matrix evolution under
   Rabi H must reproduce rabi_excited_probability()'s exact analytic
   result - validates -i[H,\rho] commutator term and RK4 integration
   independent of dissipative machinery.
4. Trace preservation: Tr(\rho) stays 1 under RK4 evolution for both damping and
   dephasing.
5. Von Neumann entropy of general density matrices via cmatrix_eigh, cross-
   checked against closed-form 2x2 formula for pure state and maximally mixed
   state.
6. Measurement collapse: deterministic outcome selection from known diagonal
   distribution, and collapse to the correct basis projector.
*/

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/lindblad.h"
#include "../physics/qubits.h"
#include "../physics/rabi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

static double density_trace(const cmatrix_t *rho) {
  double tr = 0.0;
  for (int i = 0; i < rho->nrows; i++)
    tr += CMAT(rho, i, i).re;
  return tr;
}

// Test 1: amplitude damping, H=0, exact rho11(t) = \exp(-\gamma * t)
static int test_amplitude_damping(void) {
  int n_qubits = 1;
  double gamma = 1.0;
  double dt = 1e-3;
  int steps = 1000; // T = 1.0

  cmatrix_t *H = cmatrix_alloc(2, 2); // all zero: no unitary drive
  cmatrix_t *rho = cmatrix_alloc(2, 2);
  CMAT(rho, 1, 1) = c_real(1.0); // start in |1><1| (excited)

  cmatrix_t *L = lindblad_amplitude_damping_op(n_qubits, 0, gamma);
  cmatrix_t *ops[1] = {L};

  int fail = lindblad_evolve(rho, H, ops, 1, dt, steps) != 0;

  double T = dt * steps;
  double expected_rho11 = exp(-gamma * T);
  fail |= check_close(CMAT(rho, 1, 1).re, expected_rho11, 1e-4,
                      "amplitude damping \\rho11(T=1)");
  fail |= check_close(CMAT(rho, 0, 0).re, 1.0 - expected_rho11, 1e-4,
                      "amplitude damping \\rho00(T=1)");
  fail |=
      check_close(density_trace(rho), 1.0, 1e-6, "amplitude damping Tr(\\rho)");

  cmatrix_free(H);
  cmatrix_free(rho);
  cmatrix_free(L);

  return fail;
}

// Test 2: pure dephasing, H=0, exact \rho01(t) = \rho01(0) * \exp(-\gamma * t)
static int test_dephasing(void) {
  int n_qubits = 1;
  double gamma = 1.0;
  double dt = 1e-3;
  int steps = 1000; // T = 1.0

  cmatrix_t *H = cmatrix_alloc(2, 2);
  cmatrix_t *rho = cmatrix_alloc(2, 2);
  CMAT(rho, 0, 0) = c_real(0.5);
  CMAT(rho, 0, 1) = c_real(0.5);
  CMAT(rho, 1, 0) = c_real(0.5);
  CMAT(rho, 1, 1) = c_real(0.5); // |+><+|

  cmatrix_t *L = lindblad_dephasing_op(n_qubits, 0, gamma);
  cmatrix_t *ops[1] = {L};

  int fail = lindblad_evolve(rho, H, ops, 1, dt, steps) != 0;

  double T = dt * steps;
  double expected_rho01 = 0.5 * exp(-gamma * T);
  fail |= check_close(CMAT(rho, 0, 1).re, expected_rho01, 1e-4,
                      "dephasing \\rho01(T=1)");
  fail |= check_close(CMAT(rho, 0, 0).re, 0.5, 1e-6,
                      "dephasing \\rho00 (population unaffected)");
  fail |= check_close(CMAT(rho, 1, 1).re, 0.5, 1e-6,
                      "dephasing \\rho11 (population unaffected)");
  fail |= check_close(density_trace(rho), 1.0, 1e-6, "dephasing Tr(\\rho)");

  cmatrix_free(H);
  cmatrix_free(rho);
  cmatrix_free(L);

  return fail;
}

// Test 3: closed-system (n_ops=0)
static int test_closed_system_matches_rabi(void) {
  double Omega = 2.0, Delta = 1.0;
  double dt = 1e-4;
  int steps = 5000;
  double T = dt * steps;

  cmatrix_t *H = cmatrix_alloc(2, 2);
  CMAT(H, 0, 0) = c_real(0.5 * Delta);
  CMAT(H, 0, 1) = c_real(0.5 * Omega);
  CMAT(H, 1, 0) = c_real(0.5 * Omega);
  CMAT(H, 1, 1) = c_real(-0.5 * Delta);

  cmatrix_t *rho = cmatrix_alloc(2, 2);
  CMAT(rho, 0, 0) = c_real(1.0); // start in ground state

  int fail = lindblad_evolve(rho, H, NULL, 0, dt, steps) != 0;

  double expected_pe = rabi_excited_probability(T, Omega, Delta);
  fail |= check_close(CMAT(rho, 1, 1).re, expected_pe, 1e-3,
                      "closed-system P_excited(T) vs rabi_excited_probability");
  fail |= check_close(density_trace(rho), 1.0, 1e-6, "closed-system Tr(\\rho)");
  fail |= check_close(density_purity(rho), 1.0, 1e-4,
                      "closed-system purity stays 1 (unitary)");

  cmatrix_free(H);
  cmatrix_free(rho);

  return fail;
}

// Test 4: von Neumann entropy on general (eigensolver-based) density matrices
static int test_von_neumann_entropy(void) {
  int fail = 0;

  // Pure state |0>: entropy should be 0
  cmatrix_t *rho_pure = cmatrix_alloc(2, 2);
  CMAT(rho_pure, 0, 0) = c_real(1.0);
  double S_pure = density_von_neumann_entropy(rho_pure);
  fail |= check_close(S_pure, 0.0, 1e-6, "pure state entropy");
  fail |= check_close(density_purity(rho_pure), 1.0, 1e-9, "pure state purity");

  cmatrix_free(rho_pure);

  // Maximally mixed single qubit I/2: entropy should be 1 bit
  cmatrix_t *rho_mixed = cmatrix_alloc(2, 2);
  CMAT(rho_mixed, 0, 0) = c_real(0.5);
  CMAT(rho_mixed, 1, 1) = c_real(0.5);
  double S_mixed = density_von_neumann_entropy(rho_mixed);
  fail |= check_close(S_mixed, 1.0, 1e-6, "maximally mixed qubit entropy");
  fail |= check_close(density_purity(rho_mixed), 0.5, 1e-9,
                      "maximally mixed qubit purity");

  cmatrix_free(rho_mixed);

  // Bell pair reduced to a single qubit is also maximally mixed
  int n2 = 2;
  cvector_t *psi = qstate_alloc(n2);
  qstate_apply_gate1(psi, n2, 0, hadamard_gate);
  qstate_apply_cnot(psi, n2, 0, 1);
  cmatrix_t *rho_reduced = qstate_reduced_density_single(psi, n2, 0);
  double S_closed_form = von_neumann_entropy_2x2(rho_reduced);
  double S_general = density_von_neumann_entropy(rho_reduced);
  fail |= check_close(S_general, S_closed_form, 1e-6,
                      "general vs closed-form entropy (Bell reduced state)");

  cvector_free(psi);
  cmatrix_free(rho_reduced);

  return fail;
}

// Test 5: measurement collapse against known diagonal distribution
static int test_measurement_collapse(void) {
  int fail = 0;

  cmatrix_t *rho = cmatrix_alloc(2, 2);
  CMAT(rho, 0, 0) = c_real(0.3);
  CMAT(rho, 0, 1) = c_real(0.1); // coherence
  CMAT(rho, 1, 0) = c_real(0.1);
  CMAT(rho, 1, 1) = c_real(0.7);

  cmatrix_t *rho_a = cmatrix_copy(rho);
  int outcome_a =
      density_measure_computational_basis(rho_a, 0.1); // < 0.3 -> outcome 0
  printf("  measurement u=0.10 -> outcome=%d (expected 0)\n", outcome_a);
  fail |= (outcome_a != 0);
  fail |= check_close(CMAT(rho_a, 0, 0).re, 1.0, 1e-12,
                      "collapsed \\rho00 after outcome 0");
  fail |= check_close(CMAT(rho_a, 1, 1).re, 0.0, 1e-12,
                      "collapsed \\rho11 after outcome 0");
  fail |= check_close(c_abs(CMAT(rho_a, 0, 1)), 0.0, 1e-12,
                      "collapsed coherence vanishes");

  cmatrix_free(rho_a);

  cmatrix_t *rho_b = cmatrix_copy(rho);
  int outcome_b =
      density_measure_computational_basis(rho_b, 0.5); // >= 0.3 -> outcome 1
  printf("  measurement u=0.50 -> outcome=%d (expected 1)\n", outcome_b);
  fail |= (outcome_b != 1);
  fail |= check_close(CMAT(rho_b, 1, 1).re, 1.0, 1e-12,
                      "collapsed \\rho11 after outcome 1");

  cmatrix_free(rho_b);
  cmatrix_free(rho);

  return fail;
}

int main(void) {
  int failed = 0;

  printf("Amplitude damping (T1): rho11(t) = exp(-gamma*t), H=0:\n");
  failed += test_amplitude_damping();

  printf("Pure dephasing (T2): rho01(t) = rho01(0)*exp(-gamma*t), H=0:\n");
  failed += test_dephasing();

  printf("Closed system (n_ops=0) cross-check vs rabi_excited_probability:\n");
  failed += test_closed_system_matches_rabi();

  printf("Von Neumann entropy (general eigensolver-based, cross-check vs "
         "qubits.c):\n");
  failed += test_von_neumann_entropy();

  printf("Measurement collapse (computational basis):\n");
  failed += test_measurement_collapse();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
