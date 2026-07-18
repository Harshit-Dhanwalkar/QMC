/*
Hydrogen atom: radial solver and analytic wavefunctions.
*/

#include "hydrogen.h"
#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/matrix.h"
#include "../core/special/special.h"
#include "../core/vector.h"
#include "central_potential.h"
#include "potentials.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

double hydrogen_energy_level(int n) {
  if (n < 1)
    return 0.0;
  double E = -13.6 / (n * n); // in eV
  return E * E_CHARGE;        // convert to Joules
}

eigen_t *hydrogen_radial_solve(double *r, int N, int l, double hbar,
                               double mass, double e_charge, double eps0) {
  if (!r || N < 3)
    return NULL;

  // H = -\hbar^2/(2m) d^2/dr^2 + \hbar^2/(2m) l(l+1)/r^2 - \exp^2/(4 \pi \eps0
  // r)
  double e2 = e_charge * e_charge / (4.0 * M_PI * eps0);
  return central_potential_radial_solve(r, N, l, hbar, mass, V_coulomb, &e2);
}

cvector_t *hydrogen_radial_wavefunction(double *r, int N, int n, int l) {
  if (!r || N < 1 || n < 1 || l < 0 || l >= n)
    return NULL;

  // Analytical: R_{nl}(r) = \sqrt((2 / (n a0))^3 * (n-l-1)!/(2n (n+l)!)) *
  // \exp(-r / (n a_0)) * (2r / (n a_0))^l * L_{n-l-1}^{2l+1}(2r / (n a_0))
  double a0 = AU_LENGTH;
  cvector_t *psi = cvector_alloc(N);
  if (!psi)
    return NULL;

  double norm = sqrt(pow(2.0 / (n * a0), 3) * factorial(n - l - 1) /
                     (2.0 * n * factorial(n + l)));
  for (int i = 0; i < N; i++) {
    double rho_i = 2.0 * r[i] / (n * a0);
    double L = laguerre(n - l - 1, 2 * l + 1, rho_i);
    double R = norm * exp(-rho_i / 2.0) * pow(rho_i, l) * L;
    psi->data[i] = c_real(R);
  }
  return psi;
}
