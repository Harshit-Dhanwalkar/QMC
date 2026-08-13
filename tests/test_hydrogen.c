#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../physics/hydrogen.h"
#include <stdio.h>
#include <stdlib.h>

// TODO: verify <r^2> and radial node count for n=2, l=0

int main() {
  printf(" > Testing radial wavefunction for n=1, l=0...\n");

  int N = 100;
  double r_min = 0.0, r_max = 10.0 * AU_LENGTH;
  double *r = linspace(r_min, r_max, N);
  if (!r) {
    return 1;
  }

  cvector_t *R = hydrogen_radial_wavefunction(r, N, 1, 0);
  if (!R) {
    free(r);

    return 1;
  }

  // Check normalization (\int |R|^2 * r^2 dr should be 1)
  double dr = (r_max - r_min) / (N - 1);
  double norm_sq = 0.0;
  for (int i = 0; i < N; i++) {
    double r_i = r[i];

    norm_sq += c_abs2(R->data[i]) * r_i * r_i * dr;
  }
  printf("   Normalization of n=1,l=0 radial function: %f (should be 1)\n",
         norm_sq);

  // Also test energy level: 1s = -13.6 eV
  double E1 = hydrogen_energy_level(1) / E_CHARGE;
  printf("   Hydrogen 1s energy: %f eV (expected -13.6)\n", E1);

  cvector_free(R);
  free(r);

  return 0;
}
