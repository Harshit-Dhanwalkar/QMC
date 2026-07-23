/*
 * 1D Dirac Equation (Free Particle)
 *
 * Solves the 1D Dirac equation via the complex-Hermitian real-embedding solver
 * (cmatrix_eigh_complex). For a free particle, spectrum should split into two
 * branches separated by ~2mc^2, hard-wall-grid qualitative analog of
 * theoretical continuum result.
 */

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/utils.h"
#include "../core/vector.h"
#include "../export/plot.h"
#include "../physics/relativistic.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  printf(" > 1D Dirac Equation (Free Particle)\n\n");

  int N = 50;
  double dx = 0.2;
  double hbar = 1.0, c = 1.0, m = 1.0;
  double mc2 = m * c * c;

  double *x = linspace(0.0, (N - 1) * dx, N);
  double *V = calloc(N, sizeof *V); // free particle
  if (!x || !V) {
    free(x);
    free(V);

    return 1;
  }

  printf("   Grid: N=%d, dx=%.2f, m=c=\\hbar=1 (mc^2=%.2f)\n\n", N, dx, mc2);

  eigen_t *eig = dirac_1d(x, N, V, m, hbar, c);
  if (!eig) {
    fprintf(stderr, "dirac_1d failed\n");
    free(x);
    free(V);

    return 1;
  }

  int n_pos = 0, n_neg = 0, n_gap = 0;
  for (int i = 0; i < eig->n; i++) {
    double E = eig->eigenvalues[i];
    if (E >= mc2 - 1e-6) {
      n_pos++;
    } else if (E <= -mc2 + 1e-6) {
      n_neg++;
    } else
      n_gap++;
  }
  printf("   %d states >= +mc^2 (positive-energy branch)\n", n_pos);
  printf("   %d states <= -mc^2 (negative-energy branch)\n", n_neg);
  printf("   %d states in the gap (finite-grid artifact)\n\n", n_gap);

  printf("   Lowest 5 positive-energy states:\n");
  printf("   idx    E\n");
  printf("   ---  --------\n");
  int shown = 0;
  int lowest_pos_idx = -1;
  for (int i = 0; i < eig->n && shown < 5; i++) {
    if (eig->eigenvalues[i] >= mc2 - 1e-6) {
      if (lowest_pos_idx < 0) {
        lowest_pos_idx = i;
      }
      printf("   %3d  %8.4f\n", i, eig->eigenvalues[i]);
      shown++;
    }
  }
  printf("\n");

  // Plot |upper component|^2 of lowest positive-energy state.
  // eig->eigenvectors columns are 2N-dimensional spinors: rows
  // [0,N) = upper component, [N,2N) = lower component.
  if (lowest_pos_idx >= 0) {
    cvector_t *spinor =
        cvector_from_matrix_column(eig->eigenvectors, lowest_pos_idx);
    if (spinor) {
      double *upper_prob = malloc(N * sizeof *upper_prob);

      for (int i = 0; i < N; i++) {
        upper_prob[i] = c_abs2(spinor->data[i]);
      }

      plot_opts_t opts = {0};
      opts.title = "Dirac: lowest positive-energy state (upper component)";
      opts.xlabel = "x";
      opts.ylabel = "|upper|^2";
      plot_line("dirac_lowest_positive", PLOT_FORMAT_PNG, x, upper_prob, N,
                &opts);
      printf("   Saved dirac_lowest_positive.png\n");

      free(upper_prob);
      cvector_free(spinor);
    }
  }

  eigen_free(eig);
  free(x);
  free(V);

  return 0;
}
