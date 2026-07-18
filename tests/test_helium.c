/*
Test: Helium (and helium-like ion) variational ground state.

1. Numeric golden-section minimization must match exact analytic minimizer
   (Z'_opt=Z-5/16) and analytic minimum energy
2. Helium (Z=2) must reproduce standard textbook result:
   Z'_opt=1.6875, E=-2.84765625 Hartree (~-77.5 eV).
3. Variational theorem: E_variational must be >= known experimental ground-state
   energy (-2.9037 Hartree)
4. Sanity check across several Z (H-, He, Li+, Be2+) that golden-section result
   always matches analytic one.
*/

#include "../core/constants.h"
#include "../physics/helium.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double expected, double tol,
                       const char *label) {
  double err = fabs(got - expected);
  printf("  %s: got=%.8f expected=%.8f err=%.2e\n", label, got, expected, err);

  return err > tol;
}

static int test_numeric_matches_analytic(double Z) {
  double Zeff_analytic = helium_optimal_zeff_analytic(Z);
  double E_analytic = helium_ground_state_energy_analytic(Z);

  double Zeff_numeric;
  double E_numeric =
      helium_ground_state_energy_numeric(Z, 1e-12, &Zeff_numeric);

  char label1[32], label2[32];
  snprintf(label1, sizeof label1, "Z=%.0f Z'_opt", Z);
  snprintf(label2, sizeof label2, "Z=%.0f E_opt", Z);

  // Zeff tolerance
  int fail = check_close(Zeff_numeric, Zeff_analytic, 1e-6, label1);
  fail |= check_close(E_numeric, E_analytic, 1e-8, label2);

  return fail;
}

static int test_helium_textbook_value(void) {
  double Zeff;
  double E_hartree = helium_ground_state_energy_numeric(2.0, 1e-12, &Zeff);
  double E_eV = E_hartree * (AU_ENERGY / E_CHARGE);

  printf("  Z'_opt=%.6f   (expected 1.6875)\n", Zeff);
  printf("  Z'_opt=1.6875 (theoretical values)\n");
  printf("  E = %.6f Hartree = %.4f eV\n", E_hartree, E_eV);
  printf("  E=-2.84765625 Hartree (~-77.5 eV) (theoretical values)\n");

  int fail = check_close(Zeff, 1.6875, 1e-6, "Z'_opt");
  fail |= check_close(E_hartree, -2.84765625, 1e-6, "E (Hartree)");
  fail |= check_close(E_eV, -77.5, 0.1, "E (eV)");

  return fail;
}

static int test_variational_upper_bound(void) {
  double E_variational = helium_ground_state_energy_analytic(2.0);
  double E_experimental = -2.9037; // Hartree
  printf("  E_variational=%.6f Hartree, E_experimental=%.6f Hartree\n",
         E_variational, E_experimental);
  printf("  variational theorem requires E_variational >= E_experimental\n");

  if (E_variational < E_experimental) {
    printf("  FAIL: variational energy is below true ground state, violates "
           "variational theorem\n");
    return 1;
  }
  double gap = E_variational - E_experimental;

  printf("  gap (correlation energy neglected by this ansatz): %.6f "
         "Hartree\n",
         gap);

  return (gap > 0.2) ? 1 : 0; // known gap is ~0.056 Hartree
}

int main(void) {
  int failed = 0;

  printf("Numeric golden-section vs analytic closed form (several Z):\n");
  double Z_values[4] = {1.0, 2.0, 3.0, 4.0}; // H-, He, Li+, Be2+
  for (int i = 0; i < 4; i++)
    failed += test_numeric_matches_analytic(Z_values[i]);

  printf("Helium (Z=2) theoretical value:\n");
  failed += test_helium_textbook_value();

  printf("Variational upper-bound theorem:\n");
  failed += test_variational_upper_bound();

  if (failed) {
    printf("FAILED (%d)\n", failed);
    return 1;
  }
  printf("PASS\n");
  return 0;
}
