/*
Test + demonstration: Turing-trivial vs Turing-universal for QM simulator

1. Single-qubit gates alone: Never produce entanglement, no matter how
   much superposition exists per-qubit. Demonstrated by checking
   reduced-density-matrix entropy stays exactly 0 for every qubit even
   after applying Hadamard to all of them independently. This is
   "trivial" case - state stays a product state, cost stays linear in
   n_qubits (n independent 2-dim problems), not exponential.
2. One entangling gate (CNOT) is enough to leave that regime: Bell state (H +
   CNOT) and GHZ state (H + CNOT + CNOT) checked against their exact known
   amplitudes (1/sqrt2 each on two basis states, 0 elsewhere) and their reduced
   single-qubit entropy is exactly 1 bit (maximally entangled), unambiguous
   contrast with test 1.
3. Cost is exponential by construction : state vector dimension vs n_qubits,
   printed for n=2..14.
*/

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/qubits.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

// Test 1: product state
static int test_product_state_no_entanglement(void) {
  int n_qubits = 4;
  cvector_t *psi = qstate_alloc(n_qubits);

  for (int q = 0; q < n_qubits; q++)
    qstate_apply_gate1(psi, n_qubits, q, hadamard_gate);

  int fail = 0;
  for (int q = 0; q < n_qubits; q++) {
    cmatrix_t *rho = qstate_reduced_density_single(psi, n_qubits, q);
    double S = von_neumann_entropy_2x2(rho);
    cmatrix_free(rho);
    char label[32];
    snprintf(label, sizeof label, "qubit %d entropy", q);
    fail |= check_close(S, 0.0, 1e-9, label);
  }

  cvector_free(psi);
  return fail;
}

// Test 2a: Bell state |Phi+> = (|00>+|11>)/sqrt2 via H(q0) + CNOT(0,1).
static int test_bell_state(void) {
  int n_qubits = 2;
  cvector_t *psi = qstate_alloc(n_qubits);
  qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
  qstate_apply_cnot(psi, n_qubits, 0, 1);

  double inv_sqrt2 = 1.0 / sqrt(2.0);
  int fail = 0;
  fail |= check_close(qstate_probability(psi, 0), 0.5, 1e-10, "P(|00>)");
  fail |= check_close(qstate_probability(psi, 1), 0.0, 1e-10, "P(|01>)");
  fail |= check_close(qstate_probability(psi, 2), 0.0, 1e-10, "P(|10>)");
  fail |= check_close(qstate_probability(psi, 3), 0.5, 1e-10, "P(|11>)");
  fail |= check_close(psi->data[0].re, inv_sqrt2, 1e-10, "amp(|00>).re");
  fail |= check_close(psi->data[3].re, inv_sqrt2, 1e-10, "amp(|11>).re");

  cmatrix_t *rho = qstate_reduced_density_single(psi, n_qubits, 0);
  double S = von_neumann_entropy_2x2(rho);
  cmatrix_free(rho);
  fail |= check_close(S, 1.0, 1e-9, "qubit 0 entanglement entropy (bits)");

  cvector_free(psi);
  return fail;
}

// Test 2b: GHZ state (|000>+|111>)/sqrt2, n=3, via H(q0)+CNOT(0,1)+CNOT(0,2).
static int test_ghz_state(void) {
  int n_qubits = 3;
  cvector_t *psi = qstate_alloc(n_qubits);
  qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
  qstate_apply_cnot(psi, n_qubits, 0, 1);
  qstate_apply_cnot(psi, n_qubits, 0, 2);

  double inv_sqrt2 = 1.0 / sqrt(2.0);
  int fail = 0;
  fail |= check_close(qstate_probability(psi, 0), 0.5, 1e-10, "P(|000>)");
  fail |= check_close(qstate_probability(psi, 7), 0.5, 1e-10, "P(|111>)");

  double sum_middle = 0.0;
  for (int i = 1; i < 7; i++) {
    sum_middle += qstate_probability(psi, i);
  }

  fail |= check_close(sum_middle, 0.0, 1e-10, "P(all other states)");
  fail |= check_close(psi->data[0].re, inv_sqrt2, 1e-10, "amp(|000>).re");
  fail |= check_close(psi->data[7].re, inv_sqrt2, 1e-10, "amp(|111>).re");

  for (int q = 0; q < n_qubits; q++) {
    cmatrix_t *rho = qstate_reduced_density_single(psi, n_qubits, q);
    double S = von_neumann_entropy_2x2(rho);
    cmatrix_free(rho);
    char label[32];
    snprintf(label, sizeof label, "qubit %d entropy (GHZ)", q);
    fail |= check_close(S, 1.0, 1e-9, label);
  }

  cvector_free(psi);
  return fail;
}

// Test 3: exponential state-vector growth
// Prints dimension and memory vs n_qubits so exponential-cost
static void demo_exponential_cost(void) {
  printf("  n_qubits  dimension(2^n)  memory(complex_t, MB)\n");
  for (int n = 2; n <= 14; n += 2) {
    long long dim = 1LL << n;
    double mb = (double)(dim * sizeof(complex_t)) / (1024.0 * 1024.0);
    printf("  %8d  %14lld  %8.4f\n", n, dim, mb);
  }
}

int main(void) {
  int failed = 0;

  printf("Product state (single-qubit gates only): zero entanglement:\n");
  failed += test_product_state_no_entanglement();

  printf("Bell state (H + CNOT): known amplitudes + max entanglement:\n");
  failed += test_bell_state();

  printf("GHZ state (H + CNOT + CNOT, n=3):\n");
  failed += test_ghz_state();

  printf("Exponential state-vector cost (demonstration, not a pass/fail):\n");
  demo_exponential_cost();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
