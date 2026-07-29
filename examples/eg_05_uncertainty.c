/*
 * Heisenberg Uncertainty (Minimum-Uncertainty Gaussian State)
 *
 * The harmonic-oscillator ground state
 * \psi(x) = (m * \omega/(\pi* \hbar  * ))^(1/4) * \exp(-m * \omega * x^2/(2 *
 * \hbar))
 * is theoretical minimum-uncertainty state:
 * \Delta_x * \Delta_p = \hbar/2, and <H> = \hbar * \omega/2
 */

#include "../core/constants.h"
#include "../core/utils.h"
#include "../physics/potentials.h"
#include "../physics/uncertainty.h"
#include "../physics/wavefn.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Heisenberg Uncertainty: Minimum-Uncertainty Gaussian State\n\n");

  double hbar = 1.0, m = 1.0, omega = 1.0;
  int N = 2001;
  double x_min = -8.0, x_max = 8.0;

  wavefunction_t *wf = wavefunction_alloc(N);
  if (!wf) {
    return 1;
  }
  wf->dx = (x_max - x_min) / (N - 1);

  double norm_const = pow(m * omega / (M_PI * hbar), 0.25);
  for (int i = 0; i < N; i++) {
    double xi = x_min + i * wf->dx;
    wf->x[i] = xi;
    wf->psi->data[i].re = norm_const * exp(-m * omega * xi * xi / (2.0 * hbar));
    wf->psi->data[i].im = 0.0;
  }

  uncertainty_t u = compute_uncertainties(wf);

  printf("   <x> = %.6f  (expect 0)\n", u.mean_x);
  printf("   <p> = %.6f  (expect 0)\n", u.mean_p);
  printf("   \\Delta_x = %.6f\n", u.delta_x);
  printf("   \\Delta_p = %.6f\n", u.delta_p);
  printf(
      "   \\Delta_x * \\Delta_p = %.6f  (exact minimum: \\hbar/2 = %.6f)\n\n",
      u.product, hbar / 2.0);

  double omega_param = omega;
  double E = compute_energy_expectation(wf, V_harmonic, &omega_param, m);
  printf(
      "   <H> = %.6f  (exact ground-state energy: \\hbar * \\omega/2 = %.6f)\n",
      E, hbar * omega / 2.0);

  wavefunction_free(wf);

  return 0;
}
