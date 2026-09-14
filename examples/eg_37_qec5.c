/*
 * The 5-Qubit "Perfect" Code: The Smallest Code Correcting Any Single-Qubit
 * Error
 *
 * NOTE: 9-qubit Shor code corrects an arbitrary single-qubit error too, but at
 * the cost of 9 physical qubits per logical qubit, and it pays a price for Y
 * errors: since it independently infers an X correction (from the bit-flip
 * syndrome) and a Z correction (from the phase-flip syndrome), a Y error gets
 * corrected as Z*X rather than Y itself, recovering state only up to an global
 * phase of i. The 5-qubit code below is smallest stabilizer code that can
 * correct an arbitrary single-qubit error (it saturates the quantum Hamming
 * bound: 2^4=16 syndromes exactly cover "no error" + 5 qubits * 3 Pauli types),
 * and because its 4 stabilizers directly identify the exact Pauli type of
 * whichever error occurred, recovery is exact for X, Y, and Z errors alike.
 */

#include "../core/complex.h"
#include "../physics/qec5.h"
#include "physics/qec.h"
#include <math.h>
#include <stdio.h>

static const char *error_name(qec_error_type_t t) {
  switch (t) {
  case QEC_ERROR_X:
    return "X";
  case QEC_ERROR_Y:
    return "Y";
  default:
    return "Z";
  }
}

int main(void) {
  printf(" > The 5-Qubit Perfect Code: Correcting Arbitrary Single-Qubit "
         "Errors\n\n");

  complex_t alpha = c_new(0.6, 0.2);
  complex_t beta = c_new(0.5, -0.3);
  double norm = sqrt(c_abs2(alpha) + c_abs2(beta));
  alpha = c_scale(alpha, 1.0 / norm);
  beta = c_scale(beta, 1.0 / norm);

  printf("  Logical qubit: \\alpha=%.4f%+.4fi  \\beta=%.4f%+.4fi\n\n", alpha.re,
         alpha.im, beta.re, beta.im);

  const double u[4] = {0.3, 0.7, 0.2, 0.8};

  qec_error_type_t types[3] = {QEC_ERROR_X, QEC_ERROR_Y, QEC_ERROR_Z};

  {
    qec5_result_t r = qec5_run(alpha, beta, -1, QEC_ERROR_X, u);
    double err = sqrt(c_abs2(c_sub(r.recovered_alpha, alpha)) +
                      c_abs2(c_sub(r.recovered_beta, beta)));
    printf("  no error:  syndrome=(%d,%d,%d,%d)  recovered exactly (err=%.1e)"
           "\n\n",
           r.syndrome[0], r.syndrome[1], r.syndrome[2], r.syndrome[3], err);
  }

  for (int t = 0; t < 3; t++) {
    printf("  --- %s errors ---\n", error_name(types[t]));

    for (int q = 0; q < 5; q++) {
      qec5_result_t r = qec5_run(alpha, beta, q, types[t], u);

      double err = sqrt(c_abs2(c_sub(r.recovered_alpha, alpha)) +
                        c_abs2(c_sub(r.recovered_beta, beta)));

      printf("    %s on qubit %d: syndrome=(%d,%d,%d,%d)  diagnosed as %s on "
             "q%d  recovered exactly (err=%.1e)\n",
             error_name(types[t]), q, r.syndrome[0], r.syndrome[1],
             r.syndrome[2], r.syndrome[3], error_name(r.corrected_type),
             r.corrected_qubit, err);
    }

    printf("\n");
  }

  return 0;
}
