/*
 * Quantum Phase Estimation (QPE) on a Single-Qubit Phase Gate
 *
 * NOTE: QPE (QFT/inverse-QFT) is implemented and validated in tests/test_qpe.c
 * against Reference Nielsen & Chuang's closed-form eq. 5.34 distribution. This
 * example estimates the eigenphase \phi of U = diag(1, \exp^{2 * \pi * i *
 * \phi}) acting on its own eigenstate |1>, for a \phi that does not land
 * exactly on a multiple of 1/2^n_count; showing QPE's characteristic "spread
 * across neighboring bins, sharply peaked at nearest one" behavior, then
 * repeats with more counting qubits to show the estimate's precision improving
 * as ~1/2^n_count.
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../physics/qpe.h"
#include <math.h>
#include <stdio.h>

static void run_estimate(double phi_true, int n_count) {
  cmatrix_t *U = cmatrix_alloc(2, 2);
  CMAT(U, 0, 0) = c_real(1.0);
  CMAT(U, 0, 1) = c_zero();
  CMAT(U, 1, 0) = c_zero();
  CMAT(U, 1, 1) = c_new(cos(2.0 * M_PI * phi_true), sin(2.0 * M_PI * phi_true));

  complex_t eigenstate[2] = {c_zero(), c_real(1.0)}; // |1>, U's phi-eigenstate

  qpe_distribution_t *d = qpe_estimate_phase(n_count, 1, U, eigenstate);

  printf("  n_count=%2d (%d bins): phi_estimate=%.6f  true phi=%.6f  "
         "error=%.2e  (~1/2^n_count = %.2e)\n",
         n_count, 1 << n_count, d->phi_estimate, phi_true,
         fabs(d->phi_estimate - phi_true), 1.0 / (1 << n_count));

  /* NOTE: Show the top 3 most-probable measurement outcomes, to make "spread
   * across neighboring bins" behavior visible when \phi doesn't fall exactly on
   * a representable 1/2^n_count fraction. */
  int dim = 1 << n_count;
  for (int rank = 0; rank < 3 && rank < dim; rank++) {
    int best = -1;
    double best_p = -1.0;

    for (int k = 0; k < dim; k++) {
      if (d->probabilities[k] > best_p) {
        best_p = d->probabilities[k];
        best = k;
      }
    }

    printf("    #%d most likely: bin %d/%d -> \\phi=%.6f  P=%.4f\n", rank + 1,
           best, dim, (double)best / dim, best_p);

    d->probabilities[best] = -1.0; // remove from contention for next rank
  }

  qpe_distribution_free(d);
  cmatrix_free(U);
}

int main(void) {
  printf(" > Quantum Phase Estimation: Single-Qubit Phase Gate\n\n");

  /* \phi = 3/8 is exactly representable with 3 counting qubits (QPE should be
   * exact, zero spread). \phi = 0.3 is not exactly representable at any finite
   * n_count, so QPE's estimate only converges as n_count grows. */
  printf("Case 1: \\phi = 3/8 = 0.375, exactly representable at n_count=3\n");
  run_estimate(3.0 / 8.0, 3);

  printf("\nCase 2: \\phi = 0.3, not exactly representable, estimate improves "
         "with more counting qubits\n");
  run_estimate(0.3, 3);
  run_estimate(0.3, 6);
  run_estimate(0.3, 10);

  return 0;
}
