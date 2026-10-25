/*
 * General 3D Central Potentials - Beyond Hydrogen
 */

#include "../core/matrix.h"
#include "../core/utils.h"
#include "../physics/central_potential.h"
#include "../physics/potentials.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double well_condition(double E, double V0, double a, double m,
                             double hbar) {
  double k = sqrt(2.0 * m * (E + V0)) / hbar;
  double kappa = sqrt(-2.0 * m * E) / hbar;

  return k / tan(k * a) + kappa;
}

// Run a set of 3D central-potential radial Schrodinger examples and print
// numerical spectra alongside known reference values
//
//  - isotropic harmonic oscillator energies compared to exact E(n,l)
//  - finite spherical well bound states compared to transcendental-condition
//    estimates
//
// Each section follows: build a radial grid, solve radial eigenproblem for
// selected angular momentum channels, print low-lying eigenvalues and release
// allocated resources
int main(void) {
  printf(" > General 3D Central Potentials: Beyond Hydrogen\n");

  // 3D isotropic harmonic oscillator
  // Solve for l=0..2 and compare first 3 numerical levels to exact expression
  // E(n,l)=2n+l+3/2 in units with \hbar=m=\omega=1
  printf("   3D isotropic harmonic oscillator (\\hbar=m=\\omega=1), exact "
         "E(n,l)=2n+l+3/2:\n");
  printf("   %4s  %8s  %10s  %10s\n", "l", "n", "numeric", "exact");
  {
    int N = 200;
    double *r = linspace(0.01, 8.0, N);
    double omega = 1.0;

    for (int l = 0; l <= 2; l++) {
      // Solve radial equation for fixed l using harmonic potential callback
      eigen_t *eig =
          central_potential_radial_solve(r, N, l, 1.0, 1.0, V_harmonic, &omega);
      if (eig == NULL) {
        fprintf(stderr,
                "Error: central_potential_radial_solve failed for harmonic "
                "oscillator (l=%d).\n",
                l);

        free(r);

        return 1;
      }

      // Print lowest 3 radial quantum numbers for this l
      for (int n = 0; n < 3; n++) {
        printf("   %4d  %8d  %10.5f  %10.5f\n", l, n, eig->eigenvalues[n],
               2.0 * n + l + 1.5);
      }

      eigen_free(eig);
    }

    free(r);
  }
  printf("\n");

  // 3D Finite spherical well; report l=0 bound-state energies
  printf("   3D finite spherical well (V0=5, a=3), l=0 bound states:\n");
  {
    double a = 3.0;
    double V0 = 5.0;
    int N = 220;
    double *r = linspace(0.01, 18.0, N);
    double params[2] = {a, V0};
    eigen_t *eig = central_potential_radial_solve(r, N, 0, 1.0, 1.0,
                                                  V_finite_well, params);
    if (eig == NULL) {
      fprintf(stderr, "Error: central_potential_radial_solve failed for finite "
                      "spherical well.\n");

      free(r);

      return 1;
    }

    printf("   %8s  %10s  %10s\n", "state", "numeric", "exact (root-find)");
    int npts = 4000;
    double prevE = -V0 + 1e-6;
    double prevVal = well_condition(prevE, V0, a, 1.0, 1.0);
    int found = 0;
    for (int i = 1; i < npts && found < eig->n; i++) {
      double E = -V0 + (double)i / npts * V0;
      double val = well_condition(E, V0, a, 1.0, 1.0);
      if ((val > 0) != (prevVal > 0) && fabs(val) < 100 &&
          fabs(prevVal) < 100) {
        double lo = prevE;
        double hi = E;
        int lo_sign = well_condition(lo, V0, a, 1.0, 1.0) > 0;

        for (int b = 0; b < 100; b++) {
          double mid = 0.5 * (lo + hi);
          int mid_sign = well_condition(mid, V0, a, 1.0, 1.0) > 0;
          if (lo_sign) {
            if (mid_sign) {
              lo = mid;
            } else {
              hi = mid;
            }
          } else {
            if (!mid_sign) {
              lo = mid;
            } else {
              hi = mid;
            }
          }
        }

        printf("   %8d  %10.5f  %10.5f\n", found, eig->eigenvalues[found],
               0.5 * (lo + hi));
        found++;
      }

      prevE = E;
      prevVal = val;
    }

    eigen_free(eig);
    free(r);
  }

  return 0;
}
