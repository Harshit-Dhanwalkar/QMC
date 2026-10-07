/*
Single-qubit Bloch vector on top of Lindblad solver
*/

#include "bloch.h"
#include "../core/complex.h"
#include "lindblad.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>

enum { MAX_OPS = 2 };

static int finite3(const double v[3]) {
  return isfinite(v[0]) && isfinite(v[1]) && isfinite(v[2]);
}

// Finite and inside Bloch ball (small round-off tolerance)
static int valid_bloch_vector(const double v[3]) {
  return finite3(v) &&
         sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) <= 1.0 + 1e-9;
}

int bloch_vector_from_density(const cmatrix_t *rho, double v[3]) {
  if (!rho || !v || rho->nrows != 2 || rho->ncols != 2) {
    return -1;
  }

  complex_t r01 = CMAT(rho, 0, 1);

  v[0] = 2.0 * r01.re;
  v[1] = -2.0 * r01.im;
  v[2] = CMAT(rho, 0, 0).re - CMAT(rho, 1, 1).re;

  return 0;
}

cmatrix_t *bloch_density_from_vector(const double v[3]) {
  if (!v || !valid_bloch_vector(v)) {
    return NULL;
  }

  cmatrix_t *rho = cmatrix_alloc(2, 2);
  if (!rho) {
    return NULL;
  }

  CMAT(rho, 0, 0) = c_real(0.5 * (1.0 + v[2]));
  CMAT(rho, 1, 1) = c_real(0.5 * (1.0 - v[2]));
  CMAT(rho, 0, 1) = c_new(0.5 * v[0], -0.5 * v[1]);
  CMAT(rho, 1, 0) = c_new(0.5 * v[0], 0.5 * v[1]);

  return rho;
}

int bloch_evolve(double v[3], double omega, double delta, double gamma1,
                 double gamma_phi, double dt, int steps) {
  if (!v || !valid_bloch_vector(v) || !isfinite(omega) || !isfinite(delta) ||
      !isfinite(gamma1) || gamma1 < 0.0 || !isfinite(gamma_phi) ||
      gamma_phi < 0.0 || !isfinite(dt) || dt <= 0.0 || steps < 1) {
    return -1;
  }

  cmatrix_t *rho = bloch_density_from_vector(v);
  cmatrix_t *H = cmatrix_alloc(2, 2);
  cmatrix_t *ops[MAX_OPS] = {NULL, NULL};
  int n_ops = 0;
  int rc = -2;

  if (!rho || !H) {
    goto done;
  }

  CMAT(H, 0, 0) = c_real(0.5 * delta);
  CMAT(H, 1, 1) = c_real(-0.5 * delta);
  CMAT(H, 0, 1) = c_real(0.5 * omega);
  CMAT(H, 1, 0) = c_real(0.5 * omega);

  if (gamma1 > 0.0) {
    ops[n_ops] = lindblad_amplitude_damping_op(1, 0, gamma1);
    if (!ops[n_ops]) {
      goto done;
    }

    n_ops++;
  }

  if (gamma_phi > 0.0) {
    ops[n_ops] = lindblad_dephasing_op(1, 0, gamma_phi);
    if (!ops[n_ops]) {
      goto done;
    }

    n_ops++;
  }

  if (lindblad_evolve(rho, H, n_ops ? ops : NULL, n_ops, dt, steps) != 0) {
    rc = -1;
    goto done;
  }

  rc = bloch_vector_from_density(rho, v);

done:
  cmatrix_free(rho);
  cmatrix_free(H);
  for (int i = 0; i < MAX_OPS; i++) {
    cmatrix_free(ops[i]);
  }

  return rc;
}
