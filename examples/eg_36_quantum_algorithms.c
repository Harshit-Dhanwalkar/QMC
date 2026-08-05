/*
 * Deutsch-Jozsa and Grover's Search Algorithm
 *
 * Two canonical demonstrations of quantum speedup:
 * - Deutsch-Jozsa distinguishes constant from balanced functions in a single
 *   query (classically needs up to 2^(n-1)+1 in the worst case)
 * - Grover's search finds a marked item among N with O(\sqrt(N)) queries
 *   instead of O(N).
 */

#include "../physics/quantum_algorithms.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Deutsch-Jozsa and Grover's Search Algorithm\n\n");

  printf("   --- Deutsch-Jozsa (n_input=4 bits) ---\n");
  printf("   One oracle query distinguishes constant from balanced:\n\n");

  dj_result_t r0 = deutsch_jozsa(4, DJ_CONSTANT_0, NULL, 0);
  printf("   f(x)=0 for all x:              P(all-zero)=%.4f -> %s\n",
         r0.p_all_zero, r0.is_constant ? "CONSTANT" : "balanced");

  dj_result_t r1 = deutsch_jozsa(4, DJ_CONSTANT_1, NULL, 0);
  printf("   f(x)=1 for all x:              P(all-zero)=%.4f -> %s\n",
         r1.p_all_zero, r1.is_constant ? "CONSTANT" : "balanced");

  int p1[1] = {0};
  dj_result_t rb1 = deutsch_jozsa(4, DJ_BALANCED, p1, 1);
  printf("   f(x)=x0 (single-bit parity):   P(all-zero)=%.4f -> %s\n",
         rb1.p_all_zero, rb1.is_constant ? "constant" : "BALANCED");

  int p2[3] = {0, 1, 2};
  dj_result_t rb2 = deutsch_jozsa(4, DJ_BALANCED, p2, 3);
  printf("   f(x)=x0 xor x1 xor x2:         P(all-zero)=%.4f -> %s\n\n",
         rb2.p_all_zero, rb2.is_constant ? "constant" : "BALANCED");

  printf("   --- Grover's search (n=6 qubits, N=64 states) ---\n");
  int n_qubits = 6, target = 42;
  grover_result_t r = grover_search(n_qubits, target, -1);
  printf("   Searching for state |%d> among %d states, %d Grover "
         "iterations\n",
         target, 1 << n_qubits, r.n_iterations);
  printf("   (classically would need ~%d queries on average; Grover needs "
         "O(sqrt(N))=~%d)\n",
         (1 << n_qubits) / 2, r.n_iterations);
  printf("   P(measure target) = %.4f  (vs %.4f for a uniform guess)\n\n",
         r.p_target, 1.0 / (1 << n_qubits));

  printf("   Top 5 most probable outcomes:\n");
  int *order = malloc((1 << n_qubits) * sizeof(int));
  for (int i = 0; i < (1 << n_qubits); i++) {
    order[i] = i;
  }
  for (int i = 0; i < 5; i++) {
    int best = i;
    for (int j = i + 1; j < (1 << n_qubits); j++) {
      if (r.probabilities[order[j]] > r.probabilities[order[best]]) {
        best = j;
      }
    }

    int tmp = order[i];
    order[i] = order[best];
    order[best] = tmp;
    printf("     state %2d: P=%.4f%s\n", order[i], r.probabilities[order[i]],
           order[i] == target ? "  <- target" : "");
  }

  free(order);
  free(r.probabilities);

  return 0;
}
