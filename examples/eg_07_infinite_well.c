#include "../core/constants.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../export/plot.h"
#include "/home/harshitpd/Documents/GITHUB/QMC/core/matrix.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Example 7: Infinite Square Well (Particle in a Box)
 *
 *   V(x) = 0      for 0 < x < L
 *        = \infty otherwise (hard walls: \psi(0) = \psi(L) = 0)
 *
 * Analytical:
 *   E_n = n^2 * \pi^2 * \hbar^2 / (2 * m * L^2), n = 1, 2, 3, ...
 *   \psi_n(x) = \sqrt(2/L) \sin(n * \pi * x / L)
 */

int main(void) {
  printf(" > Infinite Square Well (Particle in a Box)\n\n");

  // Parameters (natural units: m=1, hbar=1)
  double L = 1.0;
  double m = 1.0;
  double hbar = 1.0;

  int n_interior = 199; // interior grid points (walls excluded)
  double dx = L / (n_interior + 1);

  // Full grid including walls
  int n_grid = n_interior + 2;
  double *x = linspace(0.0, L, n_grid);
  if (!x) {
    fprintf(stderr, "Memory allocation failed\n");
    return 1;
  }

  printf("     Parameters:   L = %.3f, m = %.3f, hbar = %.3f\n", L, m, hbar);
  printf("     Grid: %d interior points, dx=%.6f\n\n", n_interior, dx);

  // Analytical energies for first 5 states
  printf("   Analytical energies (units hbar^2/(2mL^2)):\n");
  printf("   n   E_n\n");
  printf("  ---  -------\n");
  for (int n = 1; n <= 5; n++) {
    double E = (n * n * M_PI * M_PI * hbar * hbar) / (2.0 * m * L * L);
    printf("   %2d  %7.4f\n", n, E);
  }
  printf("\n");

  // Build Hamiltonian on interior points.
  // V=0 inside well
  double *diag = malloc(n_interior * sizeof *diag);
  double *offdiag = malloc((n_interior - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    fprintf(stderr, "Memory allocation failed\n");
    free(x);
    free(diag);
    free(offdiag);
    return 1;
  }
  double coeff = hbar * hbar / (2.0 * m * dx * dx);

  for (int i = 0; i < n_interior; i++)
    diag[i] = 2.0 * coeff;
  for (int i = 0; i < n_interior - 1; i++)
    offdiag[i] = -coeff;

  // Diagonalize
  eigen_t *eig = tridiag_eigh(diag, offdiag, n_interior);
  if (!eig) {
    fprintf(stderr, "Eigenvalue decomposition failed\n");
    free(x);
    free(diag);
    free(offdiag);
    return 1;
  }

  printf("  Lowest 5 numerical eigenvalues (units hbar^2/(2mL^2)):\n");
  printf("   n    E_num    E_ana     error %%\n");
  printf("  ---  -------  -------   --------\n");
  for (int i = 0; i < 5 && i < eig->n; i++) {
    int n = i + 1;
    double E_num = eig->eigenvalues[i];
    double E_ana = (n * n * M_PI * M_PI * hbar * hbar) / (2.0 * m * L * L);
    double err = fabs(E_num - E_ana) / E_ana * 100.0;
    printf("   %2d  %7.4f  %7.4f  %6.3f%%\n", n, E_num, E_ana, err);
  }
  printf("\n");

  // Save potential: V=0 inside the well
  double *V = calloc(n_grid, sizeof *V);
  save_potential("infinite_well_potential.dat", x, V, n_grid);

  // Save wavefunctions and plot
  for (int i = 0; i < 4 && i < eig->n; i++) {
    cvector_t *col = cvector_from_matrix_column(eig->eigenvectors, i);
    if (!col)
      continue;

    cvector_t *psi_full = cvector_alloc(n_grid);
    psi_full->data[0].re = 0.0;
    psi_full->data[0].im = 0.0;
    psi_full->data[n_grid - 1].re = 0.0;
    psi_full->data[n_grid - 1].im = 0.0;
    for (int j = 0; j < n_interior; j++)
      psi_full->data[j + 1] = col->data[j];

    // Normalize on full grid
    double norm = 0.0;
    for (int j = 0; j < n_grid; j++)
      norm += (psi_full->data[j].re * psi_full->data[j].re +
               psi_full->data[j].im * psi_full->data[j].im) *
              dx;
    double inv = (norm > 0.0) ? 1.0 / sqrt(norm) : 1.0;

    double *y = malloc(n_grid * sizeof *y);
    for (int j = 0; j < n_grid; j++) {
      psi_full->data[j].re *= inv;
      psi_full->data[j].im *= inv;
      y[j] = psi_full->data[j].re; // real eigenfunctions
    }

    char fname[64];
    snprintf(fname, sizeof fname, "infinite_well_psi_%d.dat", i + 1);
    save_wavefunction(fname, x, psi_full, n_grid);

    char plot_name[64];
    snprintf(plot_name, sizeof plot_name, "infinite_well_psi_%d", i + 1);
    plot_opts_t opts = {0};
    opts.title = "Infinite Square Well";
    opts.xlabel = "x";
    opts.ylabel = "\\psi(x)";
    opts.width = 800;
    opts.height = 600;
    plot_line(plot_name, PLOT_FORMAT_PNG, x, y, n_grid, &opts);

    free(y);
    cvector_free(psi_full);
    cvector_free(col);
    printf("    Saved %s and plot %s.png\n", fname, plot_name);
  }

  // Save eigenvalues
  int n_save = (eig->n < 10) ? eig->n : 10;
  double *E_vals = malloc(n_save * sizeof *E_vals);
  for (int i = 0; i < n_save; i++)
    E_vals[i] = eig->eigenvalues[i];
  save_eigenvalues("infinite_well_energies.dat", E_vals, n_save);
  free(E_vals);

  // Cleanup
  free(x);
  free(V);
  free(diag);
  free(offdiag);
  eigen_free(eig);

  return 0;
}
