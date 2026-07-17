/*
General 3D central-potential radial solver.

Unlocks: 3D isotropic harmonic oscillator, finite/infinite spherical
wells, Yukawa, Morse (radial), or any user-supplied potential_fn.
*/

#include "central_potential.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include "potentials.h"
#include <stdlib.h>

eigen_t *central_potential_radial_solve(double *r, int N, int l, double hbar,
                                        double mass, potential_fn V,
                                        void *params) {
  if (!r || N < 3 || !V || mass <= 0.0)
    return NULL;

  double coeff = hbar * hbar / (2.0 * mass);
  double dr = r[1] - r[0];
  if (dr <= 0.0)
    return NULL;

  double h2 = dr * dr;
  double diag_factor = 2.0 * coeff / h2;
  double offdiag_factor = -coeff / h2;

  double *diag = malloc(N * sizeof *diag);
  double *offdiag = malloc((N - 1) * sizeof *offdiag);
  if (!diag || !offdiag) {
    free(diag);
    free(offdiag);
    return NULL;
  }

  for (int i = 0; i < N; i++) {
    double r_i = r[i];
    double centrifugal =
        (r_i > 0.0) ? coeff * l * (l + 1.0) / (r_i * r_i) : 0.0;
    diag[i] = diag_factor + V(r_i, params) + centrifugal;
  }

  for (int i = 0; i < N - 1; i++)
    offdiag[i] = offdiag_factor;

  // u(r_min) = 0, u(r_max) = 0
  double boundary_val = 1e6 * diag_factor;
  diag[0] = boundary_val;
  diag[N - 1] = boundary_val;

  eigen_t *eig = tridiag_eigh(diag, offdiag, N);
  free(diag);
  free(offdiag);
  return eig;
}
