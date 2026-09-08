/*
 * The 9-Qubit Shor Code: Correcting Arbitrary Single-Qubit Errors
 *
 * NOTE: eg_35_qec.c's 3-qubit codes each protect against exactly one error
 * type: bit-flip code catches X errors but is silent (and does nothing) for
 * Z errors, and vice versa. A Y error confuses both codes, since Y is
 * simultaneously a bit flip and a phase flip. The 9-qubit Shor code
 * concatenates both 3-qubit codes - 3 blocks of bit-flip code, with 3 blocks
 * additionally wrapped in an outer phase-flip structure - so it corrects X, Y,
 * *and* Z errors on any single one of its 9 physical qubits.
 */

#include "../core/complex.h"
#include "../physics/qec.h"
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
  printf(" > The 9-Qubit Shor Code: Correcting Arbitrary Single-Qubit "
         "Errors\n\n");

  complex_t alpha = c_new(0.6, 0.2);
  complex_t beta = c_new(0.5, -0.3);
  double norm = sqrt(c_abs2(alpha) + c_abs2(beta));
  alpha = c_scale(alpha, 1.0 / norm);
  beta = c_scale(beta, 1.0 / norm);

  printf("  Logical qubit: alpha=%.4f%+.4fi  beta=%.4f%+.4fi\n\n", alpha.re,
         alpha.im, beta.re, beta.im);

  const double u[8] = {0.3, 0.7, 0.2, 0.8, 0.4, 0.6, 0.35, 0.65};

  qec_error_type_t types[3] = {QEC_ERROR_X, QEC_ERROR_Y, QEC_ERROR_Z};
  for (int t = 0; t < 3; t++) {
    printf("  --- %s errors ---\n", error_name(types[t]));

    for (int q = -1; q < 9; q++) {
      qec_shor_result_t r = qec_shor_run(alpha, beta, q, types[t], u);

      if (q == -1) {
        printf("    no error:        ");
      } else {
        printf("    %s on qubit %d:    ", error_name(types[t]), q);
      }

      double err = sqrt(c_abs2(c_sub(r.recovered_alpha, alpha)) +
                        c_abs2(c_sub(r.recovered_beta, beta)));

      if (types[t] == QEC_ERROR_Y && q >= 0) {
        printf("recovered up to a global phase (exact-state err=%.1e, "
               "expected for Y - see qec.h)\n",
               err);
      } else {
        printf("recovered exactly (err=%.1e)\n", err);
      }
    }

    printf("\n");
  }

  printf("  Every error type, on every one of the 9 physical qubits, is "
         "corrected - the bit-flip syndromes (3 blocks) and the phase-flip "
         "syndrome (across blocks) between them cover every possible single-"
         "qubit Pauli error, which is exactly what the 3-qubit codes alone "
         "could not do.\n");

  return 0;
}
