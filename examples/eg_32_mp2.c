/*
 * MP2: Second-Order Correlation Energy Beyond Hartree-Fock
 *
 * Builds directly on the closed-shell RHF SCF: MP2 uses converged HF orbitals
 * and orbital energies (occupied and virtual) to compute a second-order
 * perturbative correction for electron correlation HF itself misses entirely.
 */

#include "../core/utils.h"
#include "../physics/hartree_fock.h"
#include "../physics/mp2.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > MP2: Second-Order Correlation Energy Beyond Hartree-Fock\n\n");

  int N = 160;
  double r_min = 1e-4, r_max = 12.0;
  double *r = linspace(r_min, r_max, N);

  printf("   Helium (Z=2, 1s^2): running RHF SCF...\n");
  hf_result_t *hf = hartree_fock_atom_s_orbitals(r, N, 2.0, 1, 0.4, 1e-8, 200);
  if (!hf) {
    printf("   HF failed.\n");
    free(r);

    return 1;
  }

  double E_exact = -2.9037;
  printf("   HF converged=%d in %d iterations: E_HF = %.6f Hartree\n\n",
         hf->converged, hf->iterations, hf->total_energy);

  printf("   Adding MP2 correction, sweeping how many (s-only) virtual\n");
  printf("   orbitals are included:\n\n");
  printf("   n_virtual   E_MP2          E_HF+MP2      gap to exact\n");

  const int nv_list[6] = {2, 5, 10, 15, 20, 30};
  for (int t = 0; t < 6; t++) {
    mp2_result_t res = mp2_correlation_energy(hf, r, N, nv_list[t]);
    printf("   %8d   %+.8f   %+.6f   %.6f\n", nv_list[t], res.e_mp2,
           res.e_total, res.e_total - E_exact);
  }

  mp2_result_t best = mp2_correlation_energy(hf, r, N, 30);
  double hf_gap = hf->total_energy - E_exact;
  double mp2_gap = best.e_total - E_exact;

  printf("\n   Comparison (Hartree):\n");
  printf("     HF alone:      %.6f  (%.1f%% of the gap to exact closed)\n",
         hf->total_energy, 0.0);
  printf("     HF + MP2:      %.6f  (%.1f%% of the gap to exact closed)\n",
         best.e_total, 100.0 * (1.0 - mp2_gap / hf_gap));
  printf("     Exact:         %.6f\n\n", E_exact);
  printf("   MP2 (s-only) recovers part of the correlation energy HF misses, "
         "but not all of it");

  hf_result_free(hf);
  free(r);

  return 0;
}
