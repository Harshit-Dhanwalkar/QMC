/*
Variational Quantum Eigensolver: hardware-efficient ansatz (qubits.c) +
coordinate-descent optimizer (reusing variational.c's golden_section_minimize).
See vqe.h for the full algorithm description and known limitations.
*/

#include "vqe.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/random.h"
#include "../core/vector.h"
#include "qubits.h"
#include "variational.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

cvector_t *vqe_prepare_ansatz(int n_qubits, int n_layers, const double *theta) {
  if (n_qubits < 1 || n_layers < 1 || !theta) {
    return NULL;
  }

  cvector_t *psi = qstate_alloc(n_qubits);
  if (!psi) {
    return NULL;
  }

  complex_t g[4];

  for (int layer = 0; layer < n_layers; layer++) {
    for (int q = 0; q < n_qubits; q++) {
      // NOLINTNEXTLINE(clang-analyzer-core.CallAndMessage)
      ry_gate(theta[layer * n_qubits + q], g);
      qstate_apply_gate1(psi, n_qubits, q, g);
    }

    for (int q = 0; q < n_qubits - 1; q++) {
      qstate_apply_cnot(psi, n_qubits, q, q + 1);
    }
  }

  return psi;
}

double vqe_expectation(const cvector_t *psi, const cmatrix_t *H) {
  if (!psi || !H || H->nrows != H->ncols || H->nrows != psi->n) {
    return 0.0;
  }

  cvector_t *Hpsi = cmatrix_mv(H, psi);
  if (!Hpsi) {
    return 0.0;
  }

  double e = cvector_expect(psi, Hpsi);
  cvector_free(Hpsi);

  return e;
}

double vqe_energy(int n_qubits, int n_layers, const double *theta,
                  const cmatrix_t *H) {
  cvector_t *psi = vqe_prepare_ansatz(n_qubits, n_layers, theta);
  if (!psi) {
    return 0.0;
  }

  double e = vqe_expectation(psi, H);
  cvector_free(psi);

  return e;
}

typedef struct {
  int n_qubits;
  int n_layers;
  double *theta; // full parameter array; only index `p` varies during a single
                 // golden_section_minimize call
  int p;
  const cmatrix_t *H;
} vqe_coord_closure_t;

static double vqe_coord_eval(double x, void *params) {
  vqe_coord_closure_t *c = (vqe_coord_closure_t *)params;
  double saved = c->theta[c->p];
  c->theta[c->p] = x;
  double e = vqe_energy(c->n_qubits, c->n_layers, c->theta, c->H);
  c->theta[c->p] =
      saved; // golden_section_minimize only reads energies along way

  return e;
}

vqe_result_t vqe_run(int n_qubits, int n_layers, const cmatrix_t *H,
                     int n_sweeps, double window, uint64_t seed) {
  vqe_result_t result = {0};

  if (n_qubits < 1 || n_layers < 1 || !H || H->nrows != H->ncols ||
      H->nrows != (1 << n_qubits) || n_sweeps < 1 || window <= 0.0) {
    return result;
  }

  int n_params = n_qubits * n_layers;
  double *theta = malloc((size_t)n_params * sizeof *theta);
  if (!theta) {
    return result;
  }

  rng_state_t rng;
  rng_seed(&rng, seed);
  for (int i = 0; i < n_params; i++) {
    theta[i] = rng_uniform_range(&rng, -M_PI, M_PI);
  }

  vqe_coord_closure_t closure = {n_qubits, n_layers, theta, 0, H};

  for (int sweep = 0; sweep < n_sweeps; sweep++) {
    for (int p = 0; p < n_params; p++) {
      closure.p = p;
      double lo = theta[p] - window;
      double hi = theta[p] + window;
      double x_opt =
          golden_section_minimize(lo, hi, vqe_coord_eval, &closure, 1e-8);
      theta[p] = x_opt;
    }
  }

  result.energy = vqe_energy(n_qubits, n_layers, theta, H);
  result.theta_opt = theta;
  result.n_params = n_params;

  return result;
}

cmatrix_t *vqe_build_tfim(int n_qubits, double J, double h) {
  if (n_qubits < 1) {
    return NULL;
  }

  int dim = 1 << n_qubits;
  cmatrix_t *H = cmatrix_alloc(dim, dim);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

  int *bits = malloc((size_t)n_qubits * sizeof *bits);
  if (!bits) {
    cmatrix_free(H);

    return NULL;
  }

  for (int i = 0; i < dim; i++) {
    for (int k = 0; k < n_qubits; k++) {
      int bitpos = n_qubits - 1 - k; // qubit 0 = leftmost/MSB convention
      bits[k] = (i >> bitpos) & 1;
    }

    // Diagonal ZZ term
    double zz_sum = 0.0;
    for (int k = 0; k + 1 < n_qubits; k++) {
      double zk = bits[k] ? -1.0 : 1.0;
      double zk1 = bits[k + 1] ? -1.0 : 1.0;
      zz_sum += zk * zk1;
    }
    CMAT(H, i, i) = c_add(CMAT(H, i, i), c_real(-J * zz_sum));

    // Off-diagonal X terms: X_k flips bit k
    for (int k = 0; k < n_qubits; k++) {
      int bitpos = n_qubits - 1 - k;
      int mask = 1 << bitpos;
      int j = i ^ mask;
      CMAT(H, i, j) = c_add(CMAT(H, i, j), c_real(-h));
    }
  }

  free(bits);

  return H;
}
