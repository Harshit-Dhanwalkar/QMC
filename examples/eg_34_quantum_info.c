/*
 * Quantum Information: Teleportation, Superdense Coding, Bell's Inequality
 */

#include "../core/complex.h"
#include "../physics/quantum_info.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  printf(" > Quantum Information: Teleportation, Superdense Coding, Bell's "
         "Inequality\n\n");

  printf("   --- Quantum teleportation ---\n");
  complex_t alpha = c_real(0.6);
  complex_t beta = c_new(0.7368487952023082, 0.31153467384692046);
  printf("   Teleporting \\alpha=%.4f  \\beta=%.4f%+.4fi  (an unknown qubit "
         "state)\n",
         alpha.re, beta.re, beta.im);
  printf("   across all 4 possible Bell-measurement outcomes:\n\n");

  // forces measurement outcome 0 or 1, respectively
  const double u_for[2] = {0.1, 0.9};

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      qi_teleport_result_t r = qi_teleport(alpha, beta, u_for[i], u_for[j]);
      printf("   m1=%d m2=%d -> Bob: \\alpha=%.4f  \\beta=%.4f%+.4fi  (matches "
             "original: %s)\n",
             r.m1, r.m2, r.bob_alpha.re, r.bob_beta.re, r.bob_beta.im,
             (fabs(r.bob_alpha.re - alpha.re) < 1e-9 &&
              fabs(r.bob_beta.re - beta.re) < 1e-9)
                 ? "yes"
                 : "no");
    }
  }
  printf("\n   Every outcome, after Bob's correction, reproduces original "
         "state exactly.\n\n");

  printf("   --- Superdense coding ---\n");
  printf(
      "   2 classical bits sent using 1 qubit + a pre-shared Bell pair:\n\n");
  for (int b1 = 0; b1 <= 1; b1++) {
    for (int b2 = 0; b2 <= 1; b2++) {
      qi_superdense_result_t r = qi_superdense(b1, b2);
      printf("   sent (%d,%d) -> decoded (%d,%d)\n", b1, b2, r.decoded_bit1,
             r.decoded_bit2);
    }
  }
  printf("\n   100%% deterministic recovery, superdense coding has no failure "
         "probability by construction.\n\n");

  printf("   --- CHSH / Bell's inequality ---\n");
  double a = 0.0, a_prime = M_PI / 2.0;
  double b = M_PI / 4.0, b_prime = 3.0 * M_PI / 4.0;
  double S = qi_chsh_S(a, a_prime, b, b_prime);
  printf("   S = E(a,b) - E(a,b') + E(a',b) + E(a',b') = %.6f\n", S);
  printf("   classical (local hidden variable) bound: |S| <= 2\n");
  printf("   quantum (Tsirelson) bound:               |S| <= 2 * \\sqrt(2) = "
         "%.6f\n\n",
         2.0 * sqrt(2.0));
  printf("   The Bell state violates the classical bound by %.4f. No local "
         "hidden-variable theory can reproduce quantum mechanics predictions",
         S - 2.0);

  return 0;
}
