/*
 * Test: state-vector qubit simulator
 *
 * 1. Product state: single-qubit gates alone never entangle - reduced
 *    single-qubit entropy stays exactly 0 for every qubit even after
 *    Hadamard on all of them
 * 2. Bell state (H + CNOT): known amplitudes 1/sqrt(2) on |00> and |11>,
 *    0 elsewhere, reduced single-qubit entropy exactly 1 bit
 * 3. GHZ state (H + CNOT + CNOT, n=3): known amplitudes 1/sqrt(2) on |000>
 *    and |111>, all other amplitudes 0, reduced single-qubit entropy
 *    exactly 1 bit for every qubit
 * 4. Full-register measurement qstate_measure(): deterministic product
 *    state, Bell-state statistics, exact GHZ collapse, invalid input
 * 5. Single-qubit measurement qstate_measure_qubit(): product state
 *    uncorrelated, Bell correlation, GHZ correlation, invalid input
 *
 * A final demonstration prints the exponential state-vector cost vs
 * n_qubits (informational, not a pass/fail check)
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/random.h"
#include "../core/vector.h"
#include "../physics/qubits.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

static void test_product_state_no_entanglement(void) {
  printf("  === Test: product state (single-qubit gates only), zero "
         "entanglement ===\n");

  int n_qubits = 4;
  cvector_t *psi = qstate_alloc(n_qubits);

  for (int q = 0; q < n_qubits; q++) {
    qstate_apply_gate1(psi, n_qubits, q, hadamard_gate);
  }

  for (int q = 0; q < n_qubits; q++) {
    cmatrix_t *rho = qstate_reduced_density_single(psi, n_qubits, q);
    double S = von_neumann_entropy_2x2(rho);

    cmatrix_free(rho);

    char label[32];
    snprintf(label, sizeof label, "qubit %d entropy", q);
    check_close(S, 0.0, 1e-9, label);
  }

  cvector_free(psi);
}

static void test_bell_state(void) {
  printf("  === Test: Bell state (H + CNOT), known amplitudes and max "
         "entanglement ===\n");

  int n_qubits = 2;
  cvector_t *psi = qstate_alloc(n_qubits);
  qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
  qstate_apply_cnot(psi, n_qubits, 0, 1);

  double inv_sqrt2 = 1.0 / sqrt(2.0);
  check_close(qstate_probability(psi, 0), 0.5, 1e-10, "P(|00>)");
  check_close(qstate_probability(psi, 1), 0.0, 1e-10, "P(|01>)");
  check_close(qstate_probability(psi, 2), 0.0, 1e-10, "P(|10>)");
  check_close(qstate_probability(psi, 3), 0.5, 1e-10, "P(|11>)");
  check_close(psi->data[0].re, inv_sqrt2, 1e-10, "amp(|00>).re");
  check_close(psi->data[3].re, inv_sqrt2, 1e-10, "amp(|11>).re");

  cmatrix_t *rho = qstate_reduced_density_single(psi, n_qubits, 0);
  double S = von_neumann_entropy_2x2(rho);

  cmatrix_free(rho);

  check_close(S, 1.0, 1e-9, "qubit 0 entanglement entropy (bits)");

  cvector_free(psi);
}

static void test_ghz_state(void) {
  printf("  === Test: GHZ state (H + CNOT + CNOT, n=3) ===\n");

  int n_qubits = 3;
  cvector_t *psi = qstate_alloc(n_qubits);
  qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
  qstate_apply_cnot(psi, n_qubits, 0, 1);
  qstate_apply_cnot(psi, n_qubits, 0, 2);

  double inv_sqrt2 = 1.0 / sqrt(2.0);
  check_close(qstate_probability(psi, 0), 0.5, 1e-10, "P(|000>)");
  check_close(qstate_probability(psi, 7), 0.5, 1e-10, "P(|111>)");

  double sum_middle = 0.0;
  for (int i = 1; i < 7; i++) {
    sum_middle += qstate_probability(psi, i);
  }

  check_close(sum_middle, 0.0, 1e-10, "P(all other states)");
  check_close(psi->data[0].re, inv_sqrt2, 1e-10, "amp(|000>).re");
  check_close(psi->data[7].re, inv_sqrt2, 1e-10, "amp(|111>).re");

  for (int q = 0; q < n_qubits; q++) {
    cmatrix_t *rho = qstate_reduced_density_single(psi, n_qubits, q);
    double S = von_neumann_entropy_2x2(rho);
    cmatrix_free(rho);
    char label[32];
    snprintf(label, sizeof label, "qubit %d entropy (GHZ)", q);
    check_close(S, 1.0, 1e-9, label);
  }

  cvector_free(psi);
}

static void test_measure_deterministic_product_state(void) {
  printf("  === Test: full-register measurement on a deterministic product "
         "state ===\n");

  int n_qubits = 3;
  const double u_values[3] = {0.0, 0.37, 0.999999};

  for (int t = 0; t < 3; t++) {
    cvector_t *psi = qstate_alloc(n_qubits);
    int outcome = qstate_measure(psi, u_values[t]);

    char label[48];
    snprintf(label, sizeof label, "u=%.6f outcome", u_values[t]);
    check_close((double)outcome, 0.0, 0.0, label);

    check_close(psi->data[0].re, 1.0, 0.0, "amp[0].re == 1");
    check_close(psi->data[0].im, 0.0, 0.0, "amp[0].im == 0");

    int others_zero = 1;
    for (int i = 1; i < psi->n; i++) {
      if (psi->data[i].re != 0.0 || psi->data[i].im != 0.0) {
        others_zero = 0;
      }
    }
    check_true(others_zero, "all other amplitudes exactly zero");

    cvector_free(psi);
  }
}

static void test_measure_bell_state_statistics(void) {
  printf("  === Test: Bell-state measurement statistics (200000 trials) ===\n");

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

  printf("  counts: [%d, %d, %d, %d]  (expect ~50/50 split on 0 and 3, 1 "
         "and 2 exactly 0)\n",
         counts[0], counts[1], counts[2], counts[3]);

  check_true(counts[1] == 0, "outcome 1 never occurs (zero amplitude)");
  check_true(counts[2] == 0, "outcome 2 never occurs (zero amplitude)");

  // Binomial std err for N=200000, p=0.5 is ~0.00112; 5-\sigma ~ 0.0056
  check_close(frac0, 0.5, 0.01, "P(outcome=0) empirical frequency");
  check_close(frac3, 0.5, 0.01, "P(outcome=3) empirical frequency");
}

static void test_measure_collapse_exact_ghz(void) {
  printf("  === Test: exact collapse on GHZ state (20 trials) ===\n");

  rng_state_t rng;
  rng_seed(&rng, 777ULL);

  int n_qubits = 3;
  int bad_outcome = 0;
  int bad_amplitudes = 0;

  for (int t = 0; t < 20; t++) {
    cvector_t *psi = qstate_alloc(n_qubits);
    qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
    qstate_apply_cnot(psi, n_qubits, 0, 1);
    qstate_apply_cnot(psi, n_qubits, 0, 2);

    int outcome = qstate_measure(psi, rng_uniform(&rng));

    if (outcome != 0 && outcome != 7) {
      bad_outcome++;
    }

    for (int i = 0; i < psi->n; i++) {
      double expected_re = (i == outcome) ? 1.0 : 0.0;
      const double tol = 1e-12;
      if (fabs(psi->data[i].re - expected_re) > tol ||
          fabs(psi->data[i].im - 0.0) > tol) {
        bad_amplitudes++;
        break;
      }
    }

    cvector_free(psi);
  }

  check_true(bad_outcome == 0, "all 20 GHZ measurements return |000> or |111>");
  check_true(bad_amplitudes == 0,
             "all 20 post-measurement states are the exact collapsed basis "
             "state");
}

static void test_measure_invalid_input(void) {
  printf("  === Test: full-register measurement invalid input ===\n");

  check_true(qstate_measure(NULL, 0.5) == -1, "NULL psi rejected");
}

static void test_measure_qubit_product_state_uncorrelated(void) {
  printf("  === Test: single-qubit measurement on a product state "
         "(uncorrelated) ===\n");

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

  check_close(f_q0, 0.5, 0.01, "P(q0=0)");
  check_close(f_q1, 0.5, 0.01, "P(q1=0) after q0 measured");
}

static void test_measure_qubit_bell_correlation(void) {
  printf("  === Test: single-qubit measurement on Bell state (perfect "
         "correlation) ===\n");

  int n_qubits = 2;
  int N = 5000;
  rng_state_t rng;
  rng_seed(&rng, 2222ULL);

  int mismatches = 0;
  int bad_amplitudes = 0;

  for (int t = 0; t < N; t++) {
    cvector_t *psi = qstate_alloc(n_qubits);
    qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
    qstate_apply_cnot(psi, n_qubits, 0, 1);

    int o0 = qstate_measure_qubit(psi, n_qubits, 0, rng_uniform(&rng));
    int o1 = qstate_measure_qubit(psi, n_qubits, 1, rng_uniform(&rng));

    if (o0 != o1) {
      mismatches++;
    }

    int expected_index = (o0 << 1) | o1; // qubit0=MSB per project convention
    for (int i = 0; i < psi->n; i++) {
      double expected_re = (i == expected_index) ? 1.0 : 0.0;
      if (fabs(psi->data[i].re - expected_re) > 1e-9 ||
          fabs(psi->data[i].im) > 1e-9) {
        bad_amplitudes++;
        break;
      }
    }

    cvector_free(psi);
  }

  printf("  mismatches between q0 and q1 outcomes: %d / %d (expect 0)\n",
         mismatches, N);
  check_true(mismatches == 0, "q0 and q1 outcomes always match");
  check_true(bad_amplitudes == 0,
             "post-measurement state is the collapsed basis state");
}

static void test_measure_qubit_ghz_correlation(void) {
  printf("  === Test: single-qubit measurement on GHZ (3-qubit) state ===\n");

  int n_qubits = 3;
  int N = 5000;
  rng_state_t rng;
  rng_seed(&rng, 3333ULL);

  int mismatches = 0;
  int bad_amplitudes = 0;

  for (int t = 0; t < N; t++) {
    cvector_t *psi = qstate_alloc(n_qubits);
    qstate_apply_gate1(psi, n_qubits, 0, hadamard_gate);
    qstate_apply_cnot(psi, n_qubits, 0, 1);
    qstate_apply_cnot(psi, n_qubits, 0, 2);

    int o0 = qstate_measure_qubit(psi, n_qubits, 0, rng_uniform(&rng));
    int o1 = qstate_measure_qubit(psi, n_qubits, 1, rng_uniform(&rng));
    int o2 = qstate_measure_qubit(psi, n_qubits, 2, rng_uniform(&rng));

    if (o0 != o1 || o1 != o2) {
      mismatches++;
    }

    int expected_index = (o0 << 2) | (o1 << 1) | o2;
    for (int i = 0; i < psi->n; i++) {
      double expected_re = (i == expected_index) ? 1.0 : 0.0;
      if (fabs(psi->data[i].re - expected_re) > 1e-9 ||
          fabs(psi->data[i].im) > 1e-9) {
        bad_amplitudes++;
        break;
      }
    }

    cvector_free(psi);
  }

  printf("  mismatches among q0/q1/q2 outcomes: %d / %d (expect 0)\n",
         mismatches, N);
  check_true(mismatches == 0, "all three qubit outcomes always match");
  check_true(bad_amplitudes == 0,
             "post-measurement state is the collapsed basis state");
}

static void test_measure_qubit_invalid_input(void) {
  printf("  === Test: single-qubit measurement invalid input ===\n");

  check_true(qstate_measure_qubit(NULL, 2, 0, 0.5) == -1, "NULL psi rejected");

  cvector_t *psi = qstate_alloc(2);
  check_true(qstate_measure_qubit(psi, 2, -1, 0.5) == -1,
             "target < 0 rejected");
  check_true(qstate_measure_qubit(psi, 2, 2, 0.5) == -1,
             "target >= n_qubits rejected");
  check_true(qstate_measure_qubit(psi, 3, 0, 0.5) == -1,
             "dimension mismatch (psi->n=4 != 2^3) rejected");
  cvector_free(psi);
}

static void demo_exponential_cost(void) {
  printf("  n_qubits  dimension(2^n)  memory(complex_t, MB)\n");
  for (int n = 2; n <= 14; n += 2) {
    long long dim = 1LL << n;
    double mb = (double)(dim * sizeof(complex_t)) / (1024.0 * 1024.0);

    printf("  %8d  %14lld  %8.4f\n", n, dim, mb);
  }
}

int main(void) {
  printf(" > State-vector qubit simulator tests:\n");

  test_product_state_no_entanglement();
  test_bell_state();
  test_ghz_state();
  test_measure_deterministic_product_state();
  test_measure_bell_state_statistics();
  test_measure_collapse_exact_ghz();
  test_measure_invalid_input();
  test_measure_qubit_product_state_uncorrelated();
  test_measure_qubit_bell_correlation();
  test_measure_qubit_ghz_correlation();
  test_measure_qubit_invalid_input();

  printf("\n  === Demonstration: exponential state-vector cost ===\n");
  demo_exponential_cost();

  if (failures) {
    printf("\n%d check(s) FAILED\n", failures);
    return 1;
  }
  printf("\nAll qubit checks passed\n");

  return 0;
}
