/*
 * 3-Qubit Quantum Error Correction: Bit-Flip and Phase-Flip Codes
 *
 * Syndrome-based error correction: two ancilla qubits extract WHICH data qubit
 * (if any) errored, without ever collapsing, encoded logical qubit's
 * superposition.
 */

#include "../core/complex.h"
#include "../physics/qec.h"
#include <stdio.h>

static void demo_code(qec_code_t code, const char *name,
                      const char *error_name) {
  complex_t alpha = c_real(0.6);
  complex_t beta = c_new(0.7368487952023082, 0.31153467384692046);

  printf("   --- %s code ---\n", name);
  printf("   Logical qubit: \\alpha=%.4f  \\beta=%.4f%+.4fi\n\n", alpha.re, beta.re,
         beta.im);

  for (int e = -1; e < 3; e++) {
    qec_result_t r = qec_run(code, alpha, beta, e, 0.3, 0.7);

    if (e == -1) {
      printf("   no error injected:      ");
    } else {
      printf("   %s on qubit %d injected: ", error_name, e);
    }
    printf("syndrome=(%d,%d) -> ", r.syndrome_s3, r.syndrome_s4);

    if (r.corrected_qubit == -1) {
      printf("no correction needed   ");
    } else {
      printf("corrected qubit %d      ", r.corrected_qubit);
    }
    printf("recovered: \\alpha=%.4f \\beta=%.4f%+.4fi\n", r.recovered_alpha.re,
           r.recovered_beta.re, r.recovered_beta.im);
  }

  printf("\n");
}

int main(void) {
  printf(" > 3-Qubit Quantum Error Correction\n\n");

  demo_code(QEC_BITFLIP, "Bit-flip", "X error");
  demo_code(QEC_PHASEFLIP, "Phase-flip", "Z error");

  printf("   Every case recovers the exact original logical qubit, regardless "
         "of which data qubit (if any) the error hit syndrome measurement "
         "identifies error without ever collapsing \\alpha/\\beta\n");

  return 0;
}
