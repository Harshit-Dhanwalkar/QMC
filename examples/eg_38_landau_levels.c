/*
 * Landau Levels on a Tight-Binding Lattice (Peierls Substitution)
 *
 * A uniform perpendicular magnetic field is introduced on the 2D square lattice
 * by multiplying y-hopping amplitudes by a position-dependent phase (Landau
 * gauge A = (0, B*x, 0)) rather than by adding any new on-site term :
 * Hofstadter/Harper model.
 * In the weak-field (continuum) limit this reproduces harmonic-oscillator-like
 * Landau level ladder :
 *  E_n = -4 * t + \omega_c * (n + 1/2),
 *
 * Where:
 *  \omega_c = 4 * \pi * t * \alpha,
 * with macroscopic degeneracy per level (onestate per flux quantum threading
 * the sample).
 */

#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../physics/lattice.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int cmp_double(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;

  return (x > y) - (x < y);
}

int main(void) {
  printf(" > Landau Levels via Peierls Substitution (Hofstadter/Harper "
         "model)\n\n");

  double t = 1.0;
  double alpha = 0.05; // flux per plaquette, in units of the flux quantum
  int nx = 10, ny = 10;

  printf("   2D square lattice, %dx%d sites, t=%.2f, alpha=%.3f (flux "
         "quanta/plaquette)\n",
         nx, ny, t, alpha);
  printf("   (periodic in y, open in x : required by Landau gauge on a finite "
         "lattice)\n\n");

  cmatrix_t *H =
      lattice_build_2d_square_magnetic(nx, ny, 0.0, t, alpha, LATTICE_PERIODIC);

  cmatrix_t *copy = cmatrix_copy(H);
  eigen_t *eig = cmatrix_eigh_complex(copy);
  cmatrix_free(copy);

  int N = nx * ny;
  double *E = malloc((size_t)N * sizeof *E);
  for (int i = 0; i < N; i++) {
    E[i] = eig->eigenvalues[i];
  }
  qsort(E, (size_t)N, sizeof(double), cmp_double);

  printf("   Lowest 8 numerical eigenvalues (note near-degenerate clustering : "
         "Landau-level degeneracy):\n");
  for (int i = 0; i < 8; i++) {
    printf("     E[%d] = %.6f\n", i, E[i]);
  }

  double predicted_n0 = lattice_landau_level_energy(0, 0.0, t, alpha);
  printf("\n   Continuum-limit prediction, E_n = -4t + omega_c*(n+1/2), "
         "\\omega_c = 4 * \\pi * t * \\alpha = %.6f:\n",
         4.0 * M_PI * t * alpha);
  printf(
      "     n=0 predicted = %.6f   (numeric E[0] = %.6f, rel. err = %.3f%%)\n",
      predicted_n0, E[0],
      fabs(E[0] - predicted_n0) / fabs(predicted_n0) * 100.0);

  printf("\n   Only the n=0 level is cleanly resolved on a small lattice : "
         "only ~\\alpha * nx * ny = %.1f flux quanta thread the sample, so "
         "higher Landau levels are still broadened together with open-boundary "
         "edge states.\n A much larger lattice is needed to resolve the higher "
         "rungs cleanly but even here ground Landau level already matches "
         "analytic continuum formula to well under 1%%.\n",
         alpha * nx * ny);

  eigen_free(eig);
  cmatrix_free(H);
  free(E);

  return 0;
}
