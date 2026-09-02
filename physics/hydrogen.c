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

static const double HYDROGEN_GROUND_EV = 13.6;
static const double TWO = 2.0;
static const double FOUR = 4.0;

double hydrogen_energy_level(int n_principal) {
  if (n_principal < 1) {
    return 0.0;
  }

  double energy = -HYDROGEN_GROUND_EV / (n_principal * n_principal); // in eV

  return energy * E_CHARGE; // convert to Joules
}

eigen_t *hydrogen_radial_solve(const double *r_grid, int num_points,
                               int l_quantum, double hbar, double mass,
                               double e_charge, double eps0) {
  if (!r_grid || num_points < 3) {
    return NULL;
  }

  // H = - \hbar^2 / (2 * m) d^2/dr^2 + \hbar^2 / (2 * m) l(l + 1) / r^2
  //     - \exp^2/(4 * \pi * \eps0 * r)
  double e_squared = e_charge * e_charge / (FOUR * M_PI * eps0);

  return central_potential_radial_solve(r_grid, num_points, l_quantum, hbar,
                                        mass, V_coulomb, &e_squared);
}

cvector_t *hydrogen_radial_wavefunction(const double *r_grid, int num_points,
                                        int n_principal, int l_quantum) {
  if (!r_grid || num_points < 1 || n_principal < 1 || l_quantum < 0 ||
      l_quantum >= n_principal) {
    return NULL;
  }

  // Analytical:
  //   R_{nl}(r) = \sqrt((2 / (n a0))^3 * (n - l - 1)! / (2n * (n + l)!)) *
  //               \exp(-r / (n a_0)) * (2r / (n a_0))^l * L_{n - l - 1}^{2l +
  //               1}(2r / (n a_0))
  double bohr_radius = AU_LENGTH;
  cvector_t *psi = cvector_alloc(num_points);
  if (!psi) {
    return NULL;
  }

  double norm = sqrt(pow(TWO / (n_principal * bohr_radius), 3) *
                     factorial(n_principal - l_quantum - 1) /
                     (TWO * n_principal * factorial(n_principal + l_quantum)));
  for (int i = 0; i < num_points; i++) {
    double rho_i = TWO * r_grid[i] / (n_principal * bohr_radius);
    double laguerre_val =
        laguerre(n_principal - l_quantum - 1, 2 * l_quantum + 1, rho_i);
    double radial_val =
        norm * exp(-rho_i / TWO) * pow(rho_i, l_quantum) * laguerre_val;

    psi->data[i] = c_real(radial_val);
  }

  return psi;
}
