/*
 * Example 6: Finite Square Well - Bound States vs Well Depth
 *
 * Sweeps V_0 from 0 to V_max and finds bound state energies numerically.
 * Shows how states appear as V_0 increases - each new bound state appears
 * when V_0 crosses a threshold, visible as a kink in the E(V_0) curves.
 *
 *   V(x) = -V_0  for |x| ≤ a,   else 0
 *   Analytic threshold for n-th bound state: V_0^n = (2n-1)^2pi^2/(8a^2)
 *   (in natural units \hbar=2m=1)
 *
 * Method: finite-difference matrix eigensolver on x \in [-L, L].
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include "../export/plot.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > Finite Square Well: Bound States vs V_0\n\n");

  // Parameters
  double a = 1.0;      // half-width of well
  double L = 8.0;      // domain half-width (L >> a)
  int N = 401;         // grid points
  double V_max = 20.0; // sweep V_0 from 0 to V_max
  int n_V = 200;       // number of V_0 values
  int n_states = 5;    // track up to 5 bound states

  double dx = 2.0 * L / (N - 1);
  double *x = linspace(-L, L, N);
  if (!x)
    return 1;

  double coeff = 0.5 / (dx * dx); // \hbar^2/2m = 0.5

  // Store energy vs V0 for each state
  double *V0_arr = malloc(n_V * sizeof *V0_arr);
  double *E_arr = malloc(n_V * n_states * sizeof *E_arr);
  int *n_bound = malloc(n_V * sizeof *n_bound);
  if (!V0_arr || !E_arr || !n_bound) {
    free(x);
    return 1;
  }

  // Analytic threshold: V0_n = (2n-1)^2 * pi^2 / (8*a^2)
  printf("   Analytic thresholds for new bound states:\n");
  for (int n = 1; n <= n_states; n++) {
    double V0_thresh = (2 * n - 1) * (2 * n - 1) * M_PI * M_PI / (8.0 * a * a);
    printf("   n=%d: V_0 > %.4f\n", n, V0_thresh);
  }
  printf("\n");
  printf("   Sweeping V_0 from 0 to %.1f...\n", V_max);

  // Sweep
  for (int v = 0; v < n_V; v++) {
    double V0 = V_max * v / (n_V - 1);
    V0_arr[v] = V0;

    // Build Hamiltonian
    cmatrix_t *H = cmatrix_alloc(N, N);
    if (!H)
      continue;

    for (int i = 0; i < N; i++) {
      double Vi = (fabs(x[i]) <= a) ? -V0 : 0.0;
      CMAT(H, i, i) = c_real(2.0 * coeff + Vi);
      if (i > 0)
        CMAT(H, i, i - 1) = c_real(-coeff);
      if (i < N - 1)
        CMAT(H, i, i + 1) = c_real(-coeff);
    }

    eigen_t *eig = cmatrix_eigh(H);
    cmatrix_free(H);

    int nb = 0;
    for (int k = 0; k < n_states && k < eig->n; k++) {
      double E = eig->eigenvalues[k];
      // Bound states have E < 0
      E_arr[v * n_states + k] = E;
      if (E < 0.0)
        nb++;
    }
    n_bound[v] = nb;
    eigen_free(eig);
  }

  // Print table at selected V0 values
  printf("\n   V_0     E_1       E_2        E_3\n");
  printf("   -----  ---------  ---------  ---------\n");
  for (int v = 0; v < n_V; v += n_V / 10) {
    printf("   %5.2f", V0_arr[v]);
    for (int k = 0; k < 3; k++) {
      double E = E_arr[v * n_states + k];
      if (E < 0)
        printf("  %9.4f", E);
      else
        printf("  %9s", "---");
    }
    printf("\n");
  }

  // Save data
  {
    char path[256];
    snprintf(path, sizeof path, "%s/finite_well_energies.dat", QMC_OUTPUT_DIR);
    FILE *f = fopen(path, "w");
    if (f) {
      fprintf(f, "# V0  E1  E2  E3  E4  E5\n");
      for (int v = 0; v < n_V; v++) {
        fprintf(f, "%.6e", V0_arr[v]);
        for (int k = 0; k < n_states; k++)
          fprintf(f, "  %.6e", E_arr[v * n_states + k]);
        fprintf(f, "\n");
      }
      fclose(f);
      printf("\n   Saved finite_well_energies.dat\n");
    }
  }

  // Plot: E vs V0 for each bound state
  {
    // Extract each state's energy as a separate array
    double **E_by_state = malloc(n_states * sizeof *E_by_state);
    for (int k = 0; k < n_states; k++) {
      E_by_state[k] = malloc(n_V * sizeof **E_by_state);
      for (int v = 0; v < n_V; v++)
        E_by_state[k][v] = E_arr[v * n_states + k];
    }

    const char *labels[] = {"n=1", "n=2", "n=3", "n=4", "n=5"};
    plot_opts_t opts = {0};
    opts.title = "Finite Well: E_n vs V_0";
    opts.xlabel = "V_0";
    opts.ylabel = "E_n";
    opts.ymin = -V_max;
    opts.ymax = 1.0;
    opts.xmin = 0;
    opts.xmax = V_max;

    plot_lines("finite_well", PLOT_FORMAT_PNG, V0_arr,
               (const double **)E_by_state, n_states, n_V, labels, &opts);
    printf("   Generated finite_well.png\n");

    for (int k = 0; k < n_states; k++)
      free(E_by_state[k]);
    free(E_by_state);
  }

  free(V0_arr);
  free(E_arr);
  free(n_bound);
  free(x);
  return 0;
}
