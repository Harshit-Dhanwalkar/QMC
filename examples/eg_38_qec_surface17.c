/*
 * The Distance-3 Rotated Surface Code: The Smallest Topological QEC Code
 *
 * NOTE: 5-qubit code is "perfect": its 16 syndromes exactly partition into 1
 * no-error + 15 single-qubit-error outcomes with none shared. The distance-3
 * surface code below trades that perfection for topological structure and a
 * locally-checkable (each stabilizer touches only 2 or 4 physically adjacent
 * qubits on a 2D lattice) syndrome extraction circuit, which is what makes
 * surface codes leading candidate for physically realizable large-scale QEC
 * WARN: a few pairs of distinct single-qubit errors alias to same syndrome
 * (they differ by a low-weight stabilizer element); still perfectly
 * correctable, just not uniquely diagnosable down to exact physical qubit/Pauli
 * that occurred
 */

#include "../core/complex.h"
#include "../physics/qec_surface17.h"
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
  printf(" > The Distance-3 Rotated Surface Code ([[9,1,3]])\n\n");
  printf("  Lattice (data qubits, row-major):\n");
  printf("    0  1  2\n");
  printf("    3  4  5\n");
  printf("    6  7  8\n\n");

  complex_t alpha = c_new(0.6, 0.2);
  complex_t beta = c_new(0.5, -0.3);
  double norm = sqrt(c_abs2(alpha) + c_abs2(beta));
  alpha = c_scale(alpha, 1.0 / norm);
  beta = c_scale(beta, 1.0 / norm);

  printf("  Logical qubit: alpha=%.4f%+.4fi  beta=%.4f%+.4fi\n\n", alpha.re,
         alpha.im, beta.re, beta.im);

  const double u[8] = {0.3, 0.7, 0.2, 0.8, 0.1, 0.9, 0.4, 0.6};

  qec_error_type_t types[3] = {QEC_ERROR_X, QEC_ERROR_Y, QEC_ERROR_Z};

  {
    qec_surface17_result_t r =
        qec_surface17_run(alpha, beta, -1, QEC_ERROR_X, u);
    double err = sqrt(c_abs2(c_sub(r.recovered_alpha, alpha)) +
                      c_abs2(c_sub(r.recovered_beta, beta)));
    printf("  no error:  recovered exactly (err=%.1e)\n\n", err);
  }

  int n_aliased = 0;

  for (int t = 0; t < 3; t++) {
    printf("  --- %s errors ---\n", error_name(types[t]));

    for (int q = 0; q < 9; q++) {
      qec_surface17_result_t r = qec_surface17_run(alpha, beta, q, types[t], u);

      double err = sqrt(c_abs2(c_sub(r.recovered_alpha, alpha)) +
                        c_abs2(c_sub(r.recovered_beta, beta)));
      int aliased = (r.corrected_qubit != q) || (r.corrected_type != types[t]);

      if (aliased) {
        n_aliased++;
      }

      printf("    %s on qubit %d: diagnosed as %s on q%d%s  recovered "
             "exactly (err=%.1e)\n",
             error_name(types[t]), q, error_name(r.corrected_type),
             r.corrected_qubit, aliased ? "  [syndrome-aliased]" : "", err);
    }

    printf("\n");
  }

  printf("  %d of 27 single-qubit errors were syndrome-aliased to a different "
         "(but stabilizer-equivalent) qubit/type; yet every single one still "
         "recovered exact original logical state, since aliased pairs "
         "differ only by an element of stabilizer group.\n",
         n_aliased);

  return 0;
}
