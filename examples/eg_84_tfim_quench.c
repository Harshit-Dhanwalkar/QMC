/*
 * Quantum quench of the transverse-field Ising chain
 *
 * NOTE: A 256-site ring is prepared in ground state at h = 0.3 (ordered phase)
 * and field is switched to h = 2.0 (paramagnet), across critical point h = J
 * Free-fermion solution gives Loschmidt rate function -\ln L(t) / N, which has
 * kinks at critical times \pi (n + 1/2) / \eps_k* (dynamical quantum phase
 * transitions), and transverse magnetisation
 */

#include "../physics/tfim_quench.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  const double J = 1.0;
  const double hi = 0.3;
  const double hf = 2.0;
  const int n = 256;

  printf(" > TFIM quench h = %.1f -> %.1f, %d sites\n\n", hi, hf, n);
  printf("  critical times t*_n:");
  for (int k = 0; k < 4; k++) {
    printf(" %.4f", tfim_quench_critical_time(J, hi, hf, k));
  }

  printf("\n\n  %6s  %14s  %10s\n", "t", "rate function", "<sx>");
  for (int i = 0; i <= 24; i++) {
    double t = 0.1 * i;

    printf("  %6.2f  %14.6f  %10.5f\n", t,
           tfim_quench_loschmidt_rate(J, hi, hf, t, n),
           tfim_quench_mx(J, hi, hf, t, n));
  }

  return 0;
}
