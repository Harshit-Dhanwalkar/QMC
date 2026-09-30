/*
 * QAOA: Quantum Approximate Optimization Algorithm on MaxCut
 *
 * NOTE: eg_31_vqe.c and eg_67_vqe_pauli_hamiltonian.c minimize the energy of a
 * quantum Hamiltonian. QAOA attacks a *classical* combinatorial problem
 * instead: the cost function is diagonal in the computational basis, and a
 * depth-p circuit of alternating "phase separator" e^{-i \gamma C} and
 * "mixer" \exp^{-i \beta \sum X} layers steers uniform superposition toward
 * good bitstrings
 *
 * Problem: MaxCut on Petersen graph (10 vertices, 15 edges, triangle-free,
 * 3-regular). Its exact maximum cut is 12. At depth 1 every edge of a
 * triangle-free 3-regular graph has same closed-form expectation, so optimum is
 * <C> = |E| (1/2 + 1/(3 sqrt 3)) = 0.6925 |E| (Farhi et al. 2014), which
 * optimizer must reproduce; deeper circuits climb toward 12
 */

#include "../physics/qaoa.h"
#include <math.h>
#include <stdio.h>

#define N_VERT 10
#define N_EDGE 15

static void print_partition(int bitstring) {
  printf("    side 0: {");
  int first = 1;
  for (int q = 0; q < N_VERT; q++) {
    if (((bitstring >> (N_VERT - 1 - q)) & 1) == 0) {
      printf("%s%d", first ? "" : ",", q);
      first = 0;
    }
  }

  printf("}   side 1: {");
  first = 1;
  for (int q = 0; q < N_VERT; q++) {
    if (((bitstring >> (N_VERT - 1 - q)) & 1) == 1) {
      printf("%s%d", first ? "" : ",", q);
      first = 0;
    }
  }

  printf("}\n");
}

int main(void) {
  printf(" > QAOA: MaxCut on the Petersen graph\n\n");

  // Outer 5-cycle, five spokes, inner pentagram
  const int u[N_EDGE] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 7, 9, 6, 8};
  const int v[N_EDGE] = {1, 2, 3, 4, 0, 5, 6, 7, 8, 9, 7, 9, 6, 8, 5};

  qaoa_problem_t *pr = qaoa_problem_maxcut(N_VERT, u, v, NULL, N_EDGE);
  if (!pr) {
    printf("  failed to build problem\n");
    return 1;
  }

  int best = 0;
  double c_max = qaoa_problem_best_cost(pr, &best);
  printf("  %d vertices, %d edges, exact max cut = %.0f (brute force over "
         "2^%d bitstrings)\n\n",
         N_VERT, N_EDGE, c_max, N_VERT);

  // Depth 1: compare the optimizer to the closed form
  double p1_exact = N_EDGE * (0.5 + 1.0 / (3.0 * sqrt(3.0)));
  qaoa_result_t r1 = qaoa_run(pr, 1, 4, 6, 0.4, 2026);
  printf("  Depth-1 optimum:\n");
  printf("    QAOA <C>           = %.8f  (gamma=%.5f, beta=%.5f)\n",
         r1.expectation, r1.gamma[0], r1.beta[0]);
  printf("    closed form        = %.8f  (|E|(1/2 + 1/(3 sqrt 3)))\n",
         p1_exact);
  printf(
      "    Wang p=1 formula   = %.8f  (evaluated at the optimized angles)\n\n",
      qaoa_maxcut_p1_expectation(pr, r1.gamma[0], r1.beta[0]));

  qaoa_result_free(&r1);

  printf("  Depth sweep (4 restarts, coordinate descent + BFGS polish):\n");
  printf("    %-3s  %-10s  %-12s  %-14s\n", "p", "<C>", "approx ratio",
         "P(optimal cut)");
  qaoa_result_t last = {0};
  for (int p = 1; p <= 4; p++) {
    qaoa_result_t r = qaoa_run(pr, p, 4, 6, 0.4, 2026);
    printf("    %-3d  %-10.5f  %-12.5f  %-14.5f\n", p, r.expectation,
           r.approx_ratio, r.optimum_prob);
    if (p == 4) {
      last = r;
    } else {
      qaoa_result_free(&r);
    }
  }

  printf("\n  Most probable bitstring at p = 4: cut value %.0f of %.0f\n",
         qaoa_problem_cost_table(pr)[last.best_bitstring], c_max);
  print_partition(last.best_bitstring);

  qaoa_result_free(&last);

  // Ising ground state: same machinery, a minimization problem
  printf("\n  Ising chain (6 spins, J=1, small fields), minimization:\n");
  const int ci[5] = {0, 1, 2, 3, 4};
  const int cj[5] = {1, 2, 3, 4, 5};
  const double J[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
  const double h[6] = {0.3, -0.2, 0.1, 0.0, -0.1, 0.2};
  qaoa_problem_t *is = qaoa_problem_ising(6, h, ci, cj, J, 5, 0.0);
  if (is) {
    double e0 = qaoa_problem_best_cost(is, NULL);
    printf("    exact ground energy = %.6f\n", e0);
    for (int p = 1; p <= 3; p += 2) {
      qaoa_result_t r = qaoa_run(is, p, 6, 6, 0.4, 99);
      printf("    p=%d  <C> = %.6f  ratio = %.4f  P(ground state) = %.4f\n", p,
             r.expectation, r.approx_ratio, r.optimum_prob);

      qaoa_result_free(&r);
    }

    qaoa_problem_free(is);
  }

  qaoa_problem_free(pr);

  return 0;
}
