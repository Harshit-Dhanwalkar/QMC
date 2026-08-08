/*
 * Hydrogen Atom (radial equation)
 *
 * Solves radial Schrödinger equation for hydrogen using finite differences.
 * Compares numerical eigenvalues with analytic E_n = -13.6 eV / n^2.
 */

#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../physics/hydrogen.h"
#include "../physics/potentials.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef QMC_OUTPUT_DIR
#define QMC_OUTPUT_DIR "output"
#endif

int main(void) {
  printf(" > Hydrogen Atom (Radial Equation)\n\n");

  // Parameters (SI units)
  double hbar = HBAR;
  double mass = M_ELECTRON;
  double e_charge = E_CHARGE;
  double eps0 = EPSILON_0;

  int l = 0;   // angular momentum (s-wave)
  int N = 501; // grid points
  // double r_min = 0.0;
  double r_min = 1e-4 * AU_LENGTH; // avoid Coulomb singularity
  double r_max = 20.0 * AU_LENGTH; // 20 Bohr radii
  double dr = (r_max - r_min) / (N - 1);

  double *r = linspace(r_min, r_max, N);
  if (!r) {
    fprintf(stderr, "Memory allocation failed\n");

    return 1;
  }

  printf("   Grid: %d points from %.3e to %.3e m, dr=%.3e m\n\n", N, r_min,
         r_max, dr);

  // Solve radial equation
  eigen_t *eig = hydrogen_radial_solve(r, N, l, hbar, mass, e_charge, eps0);
  if (!eig) {
    fprintf(stderr, "Radial solve failed\n");

    free(r);
    return 1;
  }

  // Print lowest 5 energy levels (convert J -> eV)
  printf("  Lowest 5 s-states (l=0):\n");
  printf("   n   Numerical E (eV)   Analytical E (eV)   Error (%%)\n");
  printf("  ---  -----------------  ------------------  ----------\n");
  for (int i = 0; i < 5 && i < eig->n; i++) {
    double E_num = eig->eigenvalues[i] / E_CHARGE; // J -> eV
    int n = i + 1;
    double E_ana = -13.6 / (n * n);
    double err = fabs(E_num - E_ana) / fabs(E_ana) * 100.0;
    printf("   %2d   %15.6e   %15.6e   %7.2f%%\n", n, E_num, E_ana, err);
  }
  printf("\n");

  // Save radial probability densities for first 3 states
  for (int i = 0; i < 3 && i < eig->n; i++) {
    char fname[64];
    snprintf(fname, sizeof(fname), "hydrogen_radial_%d.dat", i + 1);
    cvector_t *col = cvector_from_matrix_column(eig->eigenvectors, i);

    if (col) {
      // Save r, R(r), |R(r)|^2, and probability density r^2 |R(r)|^2
      char path[512];
      snprintf(path, sizeof(path), "%s/%s", QMC_OUTPUT_DIR, fname);
      FILE *f = fopen(path, "w");
      if (f) {
        fprintf(f, "# r (m)  R(r)  |R|^2  r^2|R|^2\n");
        for (int j = 0; j < N; j++) {
          double R = col->data[j].re; // real for bound states
          double rho = R * R;

          fprintf(f, "%.6e  %.6e  %.6e  %.6e\n", r[j], R, rho,
                  r[j] * r[j] * rho);
        }

        fclose(f);
        printf("    Saved %s\n", fname);
      }

      cvector_free(col);
    }
  }

  // Cleanup
  free(r);
  eigen_free(eig);

  return 0;
}
