/*
 * General 3D Central-Potential Solver
 *
 * Demonstrates central_potential_radial_solve on 3D isotropic harmonic
 * oscillator and a finite spherical well potentials
 */

#include "../core/constants.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../export/plot.h"
#include "../physics/central_potential.h"
#include "../physics/potentials.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > General 3D Central-Potential Solver\n\n");

  // 1. 3D isotropic harmonic oscillator, l=0
  printf("   3D Isotropic Harmonic Oscillator (l=0)\n");
  printf("   Exact: E_{n_r,0} = (2 * n_r + 3/2) * \\hbar * \\omega\n\n");

  int N = 400;
  double r_min = 1e-3, r_max = 15.0;
  double *r = linspace(r_min, r_max, N);
  double hbar = 1.0, mass = 1.0, omega = 1.0;

  eigen_t *eig_ho =
      central_potential_radial_solve(r, N, 0, hbar, mass, V_harmonic, &omega);
  if (eig_ho) {
    printf("   n_r    E_num      E_exact    error\n");
    printf("   ---  --------   --------   --------\n");
    for (int k = 0; k < 4 && k < eig_ho->n; k++) {
      double E_exact = (2 * k + 1.5) * hbar * omega;
      double err = fabs(eig_ho->eigenvalues[k] - E_exact);
      printf("   %2d   %7.4f    %7.4f   %7.2e\n", k, eig_ho->eigenvalues[k],
             E_exact, err);
    }
    printf("\n");

    cvector_t *psi0 = cvector_from_matrix_column(eig_ho->eigenvectors, 0);
    if (psi0) {
      double *y = malloc(N * sizeof *y);
      for (int i = 0; i < N; i++)
        y[i] = psi0->data[i].re;

      plot_opts_t opts = {0};
      opts.title = "3D Harmonic Oscillator: ground-state u(r)";
      opts.xlabel = "r";
      opts.ylabel = "u(r)";

      // TODO: use save_wavefuntion
      plot_line("central_potential_ho_ground", PLOT_FORMAT_PNG, r, y, N, &opts);
      printf("   Saved central_potential_ho_ground.png\n\n");

      free(y);
      cvector_free(psi0);
    }

    eigen_free(eig_ho);
  }
  free(r);

  // 2. Finite spherical well
  printf("   Finite Spherical Well (l=0)\n");
  double well_params[2] = {2.0, 10.0}; // {a, V0}
  int N2 = 400;
  double *r2 = linspace(1e-3, 12.0, N2);
  eigen_t *eig_well = central_potential_radial_solve(
      r2, N2, 0, hbar, mass, V_finite_well, well_params);

  if (eig_well) {
    printf("   Bound states (E < 0):\n");
    int n_bound = 0;
    for (int k = 0; k < eig_well->n && n_bound < 5; k++) {
      if (eig_well->eigenvalues[k] < 0.0) {
        printf("   E_%d = %.4f\n", n_bound, eig_well->eigenvalues[k]);
        n_bound++;
      }
    }

    if (n_bound == 0) {
      printf("   (no bound states found for this depth/width)\n");
    }
    eigen_free(eig_well);
  }

  free(r2);

  return 0;
}
