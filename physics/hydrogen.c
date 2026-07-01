/*
Hydrogen atom: radial solver and analytic wavefunctions.
*/

#include "hydrogen.h"
#include "../core/constants.h"
#include "../core/linalg/linalg.h"
#include "../core/special/special.h"
#include "../core/matrix.h"
#include "../core/complex.h"
#include "../core/vector.h"
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
  double coeff = hbar * hbar / (2.0 * mass);
  double e2 = e_charge * e_charge / (4.0 * M_PI * eps0);

  // Build Hamiltonian: H = -coeff * d^2/dr^2 + coeff * l(l+1)/r^2 - e2/r
  // TODO: Use finite differences with logarithmic grid
  // HACK: For hydrogen, often use non-uniform grid
  double dr = r[1] - r[0];
  if (dr < 0)
    return NULL;
  double h2 = dr * dr;
  double diag_factor = 2.0 * coeff / h2;
  double offdiag_factor = -coeff / h2;

  cmatrix_t *H = cmatrix_alloc(N, N);
  if (!H)
    return NULL;

  for (int i = 0; i < N; i++) {
    double r_i = r[i];
    double V = coeff * l * (l + 1.0) / (r_i * r_i) - e2 / r_i;
    CMAT(H, i, i) = c_real(diag_factor + V);
    if (i > 0)
      CMAT(H, i, i - 1) = c_real(offdiag_factor);
    if (i < N - 1)
      CMAT(H, i, i + 1) = c_real(offdiag_factor);
  }

  // Apply boundary: R(0)=0, R(r_max)=0 (set large potential at boundaries)
  // Set diagonal at r=0 to huge value to enforce zero
  CMAT(H, 0, 0) = c_real(1e10);
  CMAT(H, N - 1, N - 1) = c_real(1e10);

  eigen_t *eig = cmatrix_eigh_generic(H);
  cmatrix_free(H);
  return eig;
}

cvector_t *hydrogen_radial_wavefunction(double *r, int N, int n, int l) {
  if (!r || N < 1 || n < 1 || l < 0 || l >= n)
    return NULL;
  // Analytical: R_{nl}(r) = sqrt((2/(n a0))^3 * (n-l-1)!/(2n (n+l)!)) *
  // exp(-r/(n a0)) * (2r/(n a0))^l * L_{n-l-1}^{2l+1}(2r/(n a0))
  double a0 = AU_LENGTH;
  cvector_t *psi = cvector_alloc(N);
  if (!psi) return NULL;

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
