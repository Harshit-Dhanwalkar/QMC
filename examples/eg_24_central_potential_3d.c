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

int main(void) {
  printf(" > General 3D Central Potentials: Beyond Hydrogen\n\n");

  // 1. 3D isotropic harmonic oscillator
  printf("   3D isotropic harmonic oscillator (\\hbar=m=\\omega=1), exact "
         "E(n,l)=2n+l+3/2:\n");
  printf("   %4s  %8s  %10s  %10s\n", "l", "n", "numeric", "exact");
  {
    int N = 200;
    double *r = linspace(0.01, 8.0, N);
    double omega = 1.0;

    for (int l = 0; l <= 2; l++) {
      eigen_t *eig =
          central_potential_radial_solve(r, N, l, 1.0, 1.0, V_harmonic, &omega);

      for (int n = 0; n < 3; n++) {
        printf("   %4d  %8d  %10.5f  %10.5f\n", l, n, eig->eigenvalues[n],
               2.0 * n + l + 1.5);
      }

      eigen_free(eig);
    }

    free(r);
  }
  printf("\n");

  // 2. 3D finite spherical well
  printf("   3D finite spherical well (V0=5, a=3), l=0 bound states:\n");
  {
    double a = 3.0, V0 = 5.0;
    int N = 220;
    double *r = linspace(0.01, 18.0, N);
    double params[2] = {a, V0};
    eigen_t *eig = central_potential_radial_solve(r, N, 0, 1.0, 1.0,
                                                  V_finite_well, params);

    printf("   %8s  %10s  %10s\n", "state", "numeric", "exact (root-find)");
    int npts = 4000;
    double prevE = -V0 + 1e-6, prevVal = well_condition(prevE, V0, a, 1.0, 1.0);
    int found = 0;
    for (int i = 1; i < npts && found < eig->n; i++) {
      double E = -V0 + (double)i / npts * V0;
      double val = well_condition(E, V0, a, 1.0, 1.0);
      if ((val > 0) != (prevVal > 0) && fabs(val) < 100 &&
          fabs(prevVal) < 100) {
        double lo = prevE, hi = E;
        double lo_sign = well_condition(lo, V0, a, 1.0, 1.0) > 0;

        for (int b = 0; b < 100; b++) {
          double mid = 0.5 * (lo + hi);
          if ((well_condition(mid, V0, a, 1.0, 1.0) > 0) == lo_sign) {
            lo = mid;
          } else {
            hi = mid;
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
