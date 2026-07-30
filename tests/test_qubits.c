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
#include "../core/random.h"
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

// Test 4a: deterministic case (product state |000>) : P(0)=1, any valid u must
// always return outcome 0, and post-measurement state must be exactly same
// basis state (nothing to collapse away)
static int test_measure_deterministic_product_state(void) {
  int n_qubits = 3;
  double u_values[3] = {0.0, 0.37, 0.999999};
  int fail = 0;

  for (int t = 0; t < 3; t++) {
    cvector_t *psi = qstate_alloc(n_qubits);
    int outcome = qstate_measure(psi, u_values[t]);

    char label[48];
    snprintf(label, sizeof label, "u=%.6f outcome", u_values[t]);
    fail |= check_close((double)outcome, 0.0, 0.0, label);

    // Post-measurement state must be exactly |000>: amp[0]=1, rest=0
    fail |= (psi->data[0].re != 1.0 || psi->data[0].im != 0.0);
    for (int i = 1; i < psi->n; i++) {
      if (psi->data[i].re != 0.0 || psi->data[i].im != 0.0) {
        fail = 1;
      }
    }

    cvector_free(psi);
  }

  return fail;
}

// Test 4b: statistical check against a known entangled state (Bell state).
// Outcomes 1 and 2 have exactly zero amplitude and must NEVER occur; 0 and
// 3 should split ~50/50 over many independent trials (fresh state rebuilt
// each trial, since qstate_measure collapses in place).
static int test_measure_bell_state_statistics(void) {
  int n_qubits = 2;
  int N = 200000;
  int counts[4] = {0, 0, 0, 0};

  rng_state_t rng;
  rng_seed(&rng, 20260730ULL);

  for (int t = 0; t < N; t++) {
    cvector_t *psi = qstate_alloc(n_qubits);
    qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
    qstate_apply_cnot(psi, n_qubits, 0, 1);

    double u = rng_uniform(&rng);
    int outcome = qstate_measure(psi, u);
    counts[outcome]++;

    cvector_free(psi);
  }

  double frac0 = (double)counts[0] / N;
  double frac3 = (double)counts[3] / N;

  printf("  counts: [%d, %d, %d, %d]  (expect ~50/50 split on 0 and 3, "
         "1 and 2 exactly 0)\n",
         counts[0], counts[1], counts[2], counts[3]);

  int fail = 0;
  fail |= (counts[1] != 0);
  fail |= (counts[2] != 0);
  // Binomial std err for N=200000, p=0.5 is ~0.00112; 5-sigma ~ 0.0056
  fail |= check_close(frac0, 0.5, 0.01, "P(outcome=0) empirical frequency");
  fail |= check_close(frac3, 0.5, 0.01, "P(outcome=3) empirical frequency");

  return fail;
}

// Test 4c: collapse correctness on a 3-outcome-capable state (GHZ, n=3) --
// after measuring, exactly one amplitude must be 1.0 and all others exactly
// 0.0, and the returned outcome must be one of the two states with nonzero
// amplitude in the original GHZ superposition (0 or 7).
static int test_measure_collapse_exact_ghz(void) {
  rng_state_t rng;
  rng_seed(&rng, 777ULL);

  int fail = 0;
  for (int t = 0; t < 20; t++) {
    int n_qubits = 3;
    cvector_t *psi = qstate_alloc(n_qubits);
    qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
    qstate_apply_cnot(psi, n_qubits, 0, 1);
    qstate_apply_cnot(psi, n_qubits, 0, 2);

    int outcome = qstate_measure(psi, rng_uniform(&rng));

    if (outcome != 0 && outcome != 7) {
      printf("  FAIL: GHZ measurement returned impossible outcome %d\n",
             outcome);
      fail = 1;
    }

    for (int i = 0; i < psi->n; i++) {
      double expected_re = (i == outcome) ? 1.0 : 0.0;
      if (psi->data[i].re != expected_re || psi->data[i].im != 0.0) {
        printf("  FAIL: post-measurement amplitude[%d] not exact\n", i);
        fail = 1;
      }
    }

    cvector_free(psi);
  }

  return fail;
}

// Test 4d: invalid input handling
static int test_measure_invalid_input(void) {
  int fail = 0;
  fail |= (qstate_measure(NULL, 0.5) != -1);

  return fail;
}

// Test 5a: single-qubit measurement on a product state (H on both qubits,
// no entanglement). Measuring qubit 0 should give ~50/50, and qubit 1's
// marginal should remain ~50/50 AFTER qubit 0 is measured, since there's
// no correlation between them in a product state -- contrast with the
// Bell-state case below.
static int test_measure_qubit_product_state_uncorrelated(void) {
  int n_qubits = 2;
  int N = 100000;
  int counts_q0[2] = {0, 0};
  int counts_q1_after[2] = {0, 0};

  rng_state_t rng;
  rng_seed(&rng, 1111ULL);

  for (int t = 0; t < N; t++) {
    cvector_t *psi = qstate_alloc(n_qubits);
    qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
    qstate_apply_gate1(psi, n_qubits, 1, hadamard_gate);

    int o0 = qstate_measure_qubit(psi, n_qubits, 0, rng_uniform(&rng));
    int o1 = qstate_measure_qubit(psi, n_qubits, 1, rng_uniform(&rng));
    counts_q0[o0]++;
    counts_q1_after[o1]++;

    cvector_free(psi);
  }

  double f_q0 = (double)counts_q0[0] / N;
  double f_q1 = (double)counts_q1_after[0] / N;
  printf("  P(q0=0)=%.4f  P(q1=0 | after q0 measured)=%.4f  (both ~0.5, "
         "uncorrelated)\n",
         f_q0, f_q1);

  int fail = 0;
  fail |= check_close(f_q0, 0.5, 0.01, "P(q0=0)");
  fail |= check_close(f_q1, 0.5, 0.01, "P(q1=0) after q0 measured");

  return fail;
}

// Test 5b: single-qubit measurement on a Bell state. Measuring qubit 0
// forces qubit 1's subsequent measurement to match EXACTLY, every time --
// this is the physically meaningful signature of entanglement, and the
// main thing distinguishing this from the product-state test above.
static int test_measure_qubit_bell_correlation(void) {
  int N = 5000;
  rng_state_t rng;
  rng_seed(&rng, 2222ULL);

  int fail = 0;
  int mismatches = 0;

  for (int t = 0; t < N; t++) {
    int n_qubits = 2;
    cvector_t *psi = qstate_alloc(n_qubits);
    qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
    qstate_apply_cnot(psi, n_qubits, 0, 1);

    int o0 = qstate_measure_qubit(psi, n_qubits, 0, rng_uniform(&rng));
    int o1 = qstate_measure_qubit(psi, n_qubits, 1, rng_uniform(&rng));

    if (o0 != o1) {
      mismatches++;
    }

    // After both qubits are measured, the state must be exactly the
    // corresponding basis state (norm 1, both amplitudes accounted for).
    int expected_index = (o0 << 1) | o1; // qubit0=MSB per project convention
    for (int i = 0; i < psi->n; i++) {
      double expected_re = (i == expected_index) ? 1.0 : 0.0;
      if (fabs(psi->data[i].re - expected_re) > 1e-9 ||
          fabs(psi->data[i].im) > 1e-9) {
        fail = 1;
      }
    }

    cvector_free(psi);
  }

  printf("  mismatches between q0 and q1 outcomes: %d / %d (expect 0)\n",
         mismatches, N);
  fail |= (mismatches != 0);

  return fail;
}

// Test 5c: invalid input handling
static int test_measure_qubit_invalid_input(void) {
  int fail = 0;
  fail |= (qstate_measure_qubit(NULL, 2, 0, 0.5) != -1);

  cvector_t *psi = qstate_alloc(2);
  fail |= (qstate_measure_qubit(psi, 2, -1, 0.5) != -1); // target < 0
  fail |= (qstate_measure_qubit(psi, 2, 2, 0.5) != -1);  // target >= n_qubits
  fail |= (qstate_measure_qubit(psi, 3, 0, 0.5) != -1);  // dimension mismatch
                                                         // (psi->n=4 != 2^3)
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

  printf("Measurement: deterministic product state:\n");
  failed += test_measure_deterministic_product_state();

  printf("Measurement: Bell state statistics (200000 trials):\n");
  failed += test_measure_bell_state_statistics();

  printf("Measurement: exact collapse on GHZ state:\n");
  failed += test_measure_collapse_exact_ghz();

  printf("Measurement: invalid input handling:\n");
  failed += test_measure_invalid_input();

  printf("Single-qubit measurement: product state, uncorrelated:\n");
  failed += test_measure_qubit_product_state_uncorrelated();

  printf("Single-qubit measurement: Bell state correlation:\n");
  failed += test_measure_qubit_bell_correlation();

  printf("Single-qubit measurement: invalid input handling:\n");
  failed += test_measure_qubit_invalid_input();

  printf("Exponential state-vector cost (demonstration, not a pass/fail):\n");
  demo_exponential_cost();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");

  return 0;
}
