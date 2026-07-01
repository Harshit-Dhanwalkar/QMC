#include "../core/complex.h"
#include "../core/constants.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../export/gnuplot_pipe.h" // if use GNUplot
#include "../physics/potentials.h"
// #include <gr.h>                   // if use GR
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int gnuplot_available(void) {
  int status = system("command -v gnuplot > /dev/null 2>&1");
  return (status == 0);
}

/*
 * Example 2: Quantum Harmonic Oscillator
 * V(x) = 0.5 * m *  \omega^2 * x^2
 *
 * Analytical:
 *   E_n = \hbar \omega (n + 1/2)
 *   \phi_n(x) = (1/sqrt(2^n n!)) * (m \omega/\pi\hbar)^(1/4)
 *            * H_n(√(m \omega/\hbar) x) * exp(-m \omega x^2 / 2\hbar)
 */

int main(void) {
  printf(" > Quantum Harmonic Oscillator\n\n");

  // Parameters (atomic units: m=1, \hbar=1, \omega=1 for simplicity)
  double omega = 1.0;
  double m = 1.0;    // set to M_ELECTRON for physical units
  double hbar = 1.0; // set HBAR for physical

  // int n_grid = 501;
  int n_grid = 51; // Grid points // TEST: different grid sizes
  double x_min = -6.0;
  double x_max = 6.0;
  double dx = (x_max - x_min) / (n_grid - 1);

  double *x = linspace(x_min, x_max, n_grid);
  if (!x) {
    fprintf(stderr, "Memory allocation failed\n");
    return 1;
  }

  printf("     Parameters:   \\omega = %.3f, m = %.3f, \\hbar = %.3f\n", omega,
         m, hbar);
  printf("     Grid: %d points from %.3f to %.3f, dx=%.6f\n\n", n_grid, x_min,
         x_max, dx);

  // Analytical solutions for first 5 states
  printf("   Analytical energies (\\hbar \\omega units):\n");
  printf("   n   E_n    <x^2> (var)\n");
  printf("  ---  -----  -----------\n");
  for (int n = 0; n < 5; n++) {
    double E = hbar * omega * (n + 0.5);
    double x2_expect = (hbar / (2 * m * omega)) * (2 * n + 1);
    printf("   %2d  %5.3f    %6.3f\n", n, E, x2_expect);
  }
  printf("\n");

  // Numerical diagonalisation of the Hamiltonian
  // H = -\hbar^2/(2m) d^2/dx^2 + 0.5 m  \omega^2 x^2
  // Discretised: H_ij = -\hbar^2/(2m) * (1/dx^2) * (δ_{i,j-1} - 2δ_{i,j} +
  // δ_{i,j+1})
  //              + 0.5 m  \omega^2 x_i^2 δ_{i,j}

  cmatrix_t *H = cmatrix_alloc(n_grid, n_grid);
  if (!H) {
    fprintf(stderr, "Memory allocation failed\n");
    free(x);
    return 1;
  }

  double coeff = -hbar * hbar / (2 * m * dx * dx); // kinetic prefactor
  double mass_omega2 = 0.5 * m * omega * omega;

  // Zero matrix
  for (int i = 0; i < n_grid; i++) {
    for (int j = 0; j < n_grid; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

  // Set diagonal and off-diagonals
  for (int i = 0; i < n_grid; i++) {
    // Kinetic: -2*coeff on diagonal, coeff on off-diagonals
    CMAT(H, i, i) = c_real(-2.0 * coeff + mass_omega2 * x[i] * x[i]);
    if (i > 0) {
      CMAT(H, i, i - 1) = c_real(coeff);
    }
    if (i < n_grid - 1) {
      CMAT(H, i, i + 1) = c_real(coeff);
    }
  }

  // Solve eigenvalue problem
  eigen_t *eig = cmatrix_eigh(H);
  if (!eig) {
    fprintf(stderr, "Eigenvalue decomposition failed\n");
    cmatrix_free(H);
    free(x);
    return 1;
  }

  printf("  Lowest 5 numerical eigenvalues (\\hbar \\omega units):\n");
  printf("   n   E_num   E_ana   error %%\n");
  printf("  ---  -----   -----   -------\n");
  for (int i = 0; i < 5 && i < eig->n; i++) {
    double E_num = eig->eigenvalues[i];
    double E_ana = hbar * omega * (i + 0.5);
    double err = fabs(E_num - E_ana) / E_ana * 100.0;
    printf("   %2d  %5.3f   %5.3f   %5.2f%%\n", i, E_num, E_ana, err);
  }
  printf("\n");

  // Save data for plotting
  // Potential
  double *V = malloc(n_grid * sizeof(double));
  for (int i = 0; i < n_grid; i++) {
    V[i] = mass_omega2 * x[i] * x[i];
  }
  save_potential("harmonic_potential.dat", x, V, n_grid);

  // First 4 eigenfunctions
  for (int i = 0; i < 4 && i < eig->n; i++) {
    char fname[64];
    snprintf(fname, sizeof(fname), "harmonic_psi_%d.dat", i);
    cvector_t *col = cvector_from_matrix_column(eig->eigenvectors, i);
    if (col) {
      save_wavefunction(fname, x, col, n_grid);
      cvector_free(col);
    }
    printf("    Saved %s\n", fname);
  }

  // Energies
  double *E_vals = malloc(10 * sizeof(double));
  for (int i = 0; i < 10 && i < eig->n; i++)
    E_vals[i] = eig->eigenvalues[i];
  save_eigenvalues("harmonic_energies.dat", E_vals,
                   (eig->n < 10) ? eig->n : 10);
  printf("\n");
  free(E_vals);

  // Plotting using GR or gnuplot
  if (!gnuplot_available()) {
    printf("GNUplot is not available - skipping plot generation.\n");
  } else {
    gnuplot_t *gp = gnuplot_open();
    if (gp) {
      gnuplot_cmd(gp, "set terminal pngcairo size 800,600");
      gnuplot_cmd(gp, "set output 'harmonic_plots.png'");
      gnuplot_cmd(gp, "set multiplot layout 2,2 title 'Harmonic Oscillator'");
      for (int i = 0; i < 4; i++) {
        char fname[64];
        snprintf(fname, sizeof(fname), "harmonic_psi_%d.dat", i);
        gnuplot_cmd(gp, "set title 'n = %d'", i);
        gnuplot_cmd(gp, "plot '%s' using 1:2 with lines title '|\\phi|^2'",
                    fname);
      }
      gnuplot_cmd(gp, "unset multiplot");
      gnuplot_cmd(gp, "set output");
      gnuplot_close(gp);
      printf("    Generated harmonic_plots.png\n");
    } else {
      fprintf(stderr,
              "Failed to open gnuplot pipe despite availability check.\n");
    }
  }

  // Cleanup
  free(x);
  free(V);
  cmatrix_free(H);
  eigen_free(eig);

  return 0;
}
