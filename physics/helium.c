/*
Helium (and helium-like ion) ground-state energy via classic
effective-nuclear-charge variational method derivation.
*/

#include "helium.h"
#include "variational.h"
#include <math.h>

double helium_variational_energy(double Z_eff, double Z) {
  return Z_eff * Z_eff - 2.0 * Z * Z_eff + 0.625 * Z_eff;
}

double helium_optimal_zeff_analytic(double Z) { return Z - 5.0 / 16.0; }

double helium_ground_state_energy_analytic(double Z) {
  double d = Z - 5.0 / 16.0;

  return -(d * d);
}

static double helium_energy_adapter(double Z_eff, void *params) {
  double Z = *(double *)params;

  return helium_variational_energy(Z_eff, Z);
}

double helium_ground_state_energy_numeric(double Z, double tol,
                                          double *Zeff_opt_out) {
  double lo = (0.1 * Z > 0.01) ? 0.1 * Z : 0.01;
  double hi = 2.0 * Z;
  double Zeff_opt =
      golden_section_minimize(lo, hi, helium_energy_adapter, &Z, tol);

  if (Zeff_opt_out) {
    *Zeff_opt_out = Zeff_opt;
  }

  return helium_variational_energy(Zeff_opt, Z);
}
