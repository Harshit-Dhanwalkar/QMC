/*
 * Driven Two-Level Systems - Beyond Time-Independent, Beyond RWA
 *
 * NOTE: rabi.c handles a time-independent drive, in closed form, using
 * rotating-wave approximation (RWA). This module handles both : \Delta(t) and
 * \Omega(t) can be arbitrary caller-supplied functions of time (chirped pulses,
 * linear sweeps, ...), and driven_two_level_evolve_lab_frame goes back to
 * actual lab-frame Hamiltonian, keeping counter-rotating term RWA throws away.
 */

#include "../core/complex.h"
#include "../core/vector.h"
#include "../physics/driven.h"
#include "../physics/rabi.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  printf(" > Driven Two-Level Systems: Sweeps and Beyond RWA\n\n");

  // 1. Landau-Zener: diabatic vs adiabatic passage through crossing
  printf("   Landau-Zener sweep \\Delta(t)=\\alpha*t through coupling "
         "Omega=1.0:\n");
  printf("   %10s  %14s  %14s\n", "\\alpha", "P(diabatic)", "formula");
  {
    double Omega = 1.0;
    double alphas[] = {0.1, 0.3, 1.0, 3.0, 20.0};
    double Ts[] = {200.0, 40.0, 20.0, 15.0, 10.0};
    double dts[] = {2e-3, 2e-3, 5e-4, 2e-4, 1e-4};

    for (int i = 0; i < 5; i++) {
      double alpha = alphas[i], T = Ts[i], dt = dts[i];
      int steps = (int)(2.0 * T / dt);
      cvector_t *psi = cvector_alloc(2);
      psi->data[0] = c_real(1.0);
      driven_two_level_evolve(psi, time_fn_linear_ramp, &alpha,
                              time_fn_constant, &Omega, -T, dt, steps);
      double p_diabatic = c_abs2(psi->data[0]);
      printf("   %10.2f  %14.6f  %14.6e\n", alpha, p_diabatic,
             landau_zener_probability(Omega, alpha));
      cvector_free(psi);
    }
    printf("   (- Small \\alpha = slow/adiabatic sweep -> follows ground state "
           "across crossing, ends up flipped.\n");
    printf("    - Large \\alpha = fast/diabatic sweep -> no time to react, "
           "stays put)\n");
  }
  printf("\n");

  // 2. Beyond RWA: Bloch-Siegert shift grows with drive strength
  printf("   Lab-frame (no RWA) vs RWA prediction, resonant drive, growing "
         "\\Omega_0/\\omega_0 (weak-driving assumption behind RWA breaks down "
         "as this ratio grows):\n");
  printf("   %10s  %14s  %14s  %10s\n", "\\Omega0/\\omega0", "P_lab (exact)",
         "P_RWA", "|diff|");
  {
    double omega0 = 50.0, omega_L = 50.0, T = 3.0, dt = 1e-4;
    int steps = (int)(T / dt);
    double ratios[] = {0.02, 0.1, 0.3, 0.6};

    for (int i = 0; i < 4; i++) {
      double Omega0 = ratios[i] * omega0;
      cvector_t *psi = cvector_alloc(2);
      psi->data[0] = c_real(1.0);
      driven_two_level_evolve_lab_frame(psi, omega0, Omega0, omega_L, 0.0, 0.0,
                                        dt, steps);
      double p_lab = c_abs2(psi->data[1]);
      // RWA-equivalent parameters for THIS convention: Omega_rwa = Omega0
      // (no extra 1/2 -- see driven.h's docstring for why).
      double p_rwa = rabi_excited_probability(T, Omega0, omega_L - omega0);
      printf("   %10.2f  %14.6f  %14.6f  %10.2e\n", ratios[i], p_lab, p_rwa,
             fabs(p_lab - p_rwa));
      cvector_free(psi);
    }
    printf("   (RWA agrees closely while \\Omega0 << \\omega0, and visibly "
           "drifts away. The Bloch-Siegert shift as drive gets strong enough "
           "to make counter-rotating term matter)\n");
  }

  return 0;
}
