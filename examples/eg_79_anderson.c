/*
 * Anderson localisation on a 1D chain
 *
 * NOTE: A particle started on one site of a tight-binding chain (hop t = 1)
 * spreads ballistically, <x^2> = 2 t^2 tau^2, until random on-site energies
 * in [-W/2, W/2] are added. Then it stops: width saturates at a few
 * localisation lengths and inverse participation ratio stays finite
 *
 * Prints rms width and IPR against time for several disorder strengths
 */

#include "../physics/anderson.h"
#include "complex.h"
#include "vector.h"
#include <math.h>
#include <stdio.h>

enum { N = 801 };

int main(void) {
  static double eps[N];
  const double strengths[] = {0.0, 1.0, 3.0, 8.0};

  printf(" > Anderson localisation, N = %d sites, t = 1\n\n", N);
  printf("  %6s | %7s %9s %7s | %7s %9s %7s\n", "W", "\\tau", "width", "IPR",
         "\\tau", "width", "IPR");

  for (int s = 0; s < 4; s++) {
    cvector_t *psi = cvector_alloc(N);
    if (!psi) {
      return 1;
    }

    anderson_disorder(eps, N, strengths[s], 2024);
    cvector_fill(psi, c_new(0.0, 0.0));
    psi->data[N / 2] = c_new(1.0, 0.0);
    double line[3][3];

    for (int k = 0; k < 3; k++) {
      // \tau = 10, 40, 100
      int steps = (k == 0) ? 1000 : (k == 1 ? 3000 : 6000);

      if (anderson_evolve(psi, eps, 1.0, 0.01, steps) != 0) {
        cvector_free(psi);

        return 1;
      }

      line[k][0] = (k == 0) ? 10.0 : (k == 1 ? 40.0 : 100.0);
      line[k][1] = anderson_width(psi);
      line[k][2] = anderson_ipr(psi);
    }

    for (int k = 0; k < 3; k++) {
      printf("  %6.1f | %7.0f %9.3f %7.4f\n", strengths[s], line[k][0],
             line[k][1], line[k][2]);
    }

    printf("\n");

    cvector_free(psi);
  }

  printf("  Weak-disorder localisation length at E = 1: xi = %.1f sites for "
         "W = 1\n", 1.0 / anderson_lyapunov_weak(1.0, 1.0, 1.0));

  return 0;
}
