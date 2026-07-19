/*
 * Quantum Harmonic Oscillator
 * V(x) = 0.5 * m * \omega^2 * x^2
 *
 * Analytical:
 *   E_n = \hbar \omega (n + 1/2)
 *   \psi_n(x) = (1 / \sqrt(2^n n!)) * (m * \omega / \pi * \hbar)^(1/4)
 *              * H_n(\sqrt(m * \omega / \hbar) x)
 *              * \exp(-m * \omega x^2 / 2 * \hbar)
 */

#include "../core/constants.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../export/plot.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Quantum Harmonic Oscillator\n\n");

  // Parameters (atomic units: m=1, \hbar=1, \omega=1)
  double omega = 1.0;
  double m = 1.0;
  double hbar = 1.0;

  int n_grid = 51;
  double x_min = -6.0;
  double x_max = 6.0;
  double dx = (x_max - x_min) / (n_grid - 1);

  double *x = linspace(x_min, x_max, n_grid);
  if (!x) {
    fprintf(stderr, "Memory allocation failed\n");
    return 1;
  }

  printf("     Parameters:   omega = %.3f, m = %.3f, hbar = %.3f\n", omega, m,
         hbar);
  printf("     Grid: %d points from %.3f to %.3f, dx=%.6f\n\n", n_grid, x_min,
         x_max, dx);

  // Analytical energies for first 5 states
  printf("   Analytical energies (\\hbar \\omega units):\n");
  printf("   n   E_n    <x^2>\n");
  printf("  ---  -----  -----\n");
  for (int n = 0; n < 5; n++) {
    double E = hbar * omega * (n + 0.5);
    double x2_expect = (hbar / (2 * m * omega)) * (2 * n + 1);
    printf("   %2d  %5.3f   %6.3f\n", n, E, x2_expect);
  }
  printf("\n");

  // Build Hamiltonian (tridiagonal): diag[i] + offdiag connecting
  // neighbors
  double *diag = malloc(n_grid * sizeof *diag);
  double *offdiag = malloc((n_grid - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    fprintf(stderr, "Memory allocation failed\n");
    free(x);
    free(diag);
    free(offdiag);
    return 1;
  }
  double coeff = -hbar * hbar / (2 * m * dx * dx);
  double mass_omega2 = 0.5 * m * omega * omega;

  for (int i = 0; i < n_grid; i++)
    diag[i] = -2.0 * coeff + mass_omega2 * x[i] * x[i];
  for (int i = 0; i < n_grid - 1; i++)
    offdiag[i] = coeff;

  // Diagonalize
  eigen_t *eig = tridiag_eigh(diag, offdiag, n_grid);
  if (!eig) {
    fprintf(stderr, "Eigenvalue decomposition failed\n");
    free(x);
    free(diag);
    free(offdiag);
    return 1;
  }

  printf("  Lowest 5 numerical eigenvalues (hbar omega units):\n");
  printf("   n   E_num   E_ana   error %%\n");
  printf("  ---  -----   -----   -------\n");
  for (int i = 0; i < 5 && i < eig->n; i++) {
    double E_num = eig->eigenvalues[i];
    double E_ana = hbar * omega * (i + 0.5);
    double err = fabs(E_num - E_ana) / E_ana * 100.0;
    printf("   %2d  %5.3f   %5.3f   %5.2f%%\n", i, E_num, E_ana, err);
  }
  printf("\n");

  // Save potential
  double *V = malloc(n_grid * sizeof(double));
  for (int i = 0; i < n_grid; i++)
    V[i] = mass_omega2 * x[i] * x[i];
  save_potential("harmonic_potential.dat", x, V, n_grid);

  // Save wavefunctions and plot
  for (int i = 0; i < 4 && i < eig->n; i++) {
    cvector_t *col = cvector_from_matrix_column(eig->eigenvectors, i);
    if (!col)
      continue;

    double *y = malloc(n_grid * sizeof(double));
    for (int j = 0; j < n_grid; j++)
      y[j] = col->data[j].re; // real eigenfunctions

    // Save data file (real part)
    char fname[64];
    snprintf(fname, sizeof(fname), "harmonic_psi_%d.dat", i);
    save_wavefunction(fname, x, col, n_grid);

    char plot_name[64];
    snprintf(plot_name, sizeof(plot_name), "psi_%d", i);
    plot_opts_t opts = {0};
    opts.title = "Harmonic Oscillator";
    opts.xlabel = "x";
    opts.ylabel = "ψ(x)";
    opts.width = 800;
    opts.height = 600;
    plot_line(plot_name, PLOT_FORMAT_PNG, x, y, n_grid, &opts);

    free(y);
    cvector_free(col);
    printf("    Saved %s and plot %s.png\n", fname, plot_name);
  }

  // Save eigenvalues
  double *E_vals = malloc(10 * sizeof(double));
  for (int i = 0; i < 10 && i < eig->n; i++)
    E_vals[i] = eig->eigenvalues[i];
  save_eigenvalues("harmonic_energies.dat", E_vals,
                   (eig->n < 10) ? eig->n : 10);
  free(E_vals);

  // Cleanup
  free(x);
  free(V);
  free(diag);
  free(offdiag);
  eigen_free(eig);

  return 0;
}
