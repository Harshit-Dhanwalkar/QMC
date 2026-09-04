/*
Noisy VQE: drives lindblad.c's existing GKSL machinery gate-by-gate to simulate
the vqe.c hardware-efficient ansatz under T1/T2 noise, plus a
zero-noise-extrapolation (ZNE) mitigation estimator.
*/

#include "vqe_noisy.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include "lindblad.h"
#include "qubits.h"
#include <math.h>
#include <stdlib.h>

// \rho <- U \rho U^\dagger, freeing input \rho and returning a new matrix
static cmatrix_t *dm_conjugate_free(const cmatrix_t *U, cmatrix_t *rho) {
  cmatrix_t *Urho = cmatrix_multiply(U, rho);
  cmatrix_free(rho);
  if (!Urho) {
    return NULL;
  }

  cmatrix_t *Udag = cmatrix_adjoint(U);
  if (!Udag) {
    cmatrix_free(Urho);

    return NULL;
  }

  cmatrix_t *result = cmatrix_multiply(Urho, Udag);
  cmatrix_free(Urho);
  cmatrix_free(Udag);

  return result;
}

/*
 * NOTE: CNOT conjugation on a density matrix via the permutation formula:
 *  \rho'[i][j] = \rho[perm(i)][perm(j)]
 *
 * Where
 *  perm flips target's bit whenever control's bit is 1 (0-indexed, qubit 0 =
 *  leftmost/MSB).
 */
static cmatrix_t *dm_apply_cnot_free(cmatrix_t *rho, int n_qubits, int control,
                                     int target) {
  int dim = rho->nrows;
  int bitpos_c = n_qubits - 1 - control;
  int bitpos_t = n_qubits - 1 - target;

  cmatrix_t *out = cmatrix_alloc(dim, dim);
  if (!out) {
    cmatrix_free(rho);

    return NULL;
  }

  for (int i = 0; i < dim; i++) {
    int pi = ((i >> bitpos_c) & 1) ? (i ^ (1 << bitpos_t)) : i;

    for (int j = 0; j < dim; j++) {
      int pj = ((j >> bitpos_c) & 1) ? (j ^ (1 << bitpos_t)) : j;

      CMAT(out, i, j) = CMAT(rho, pi, pj);
    }
  }

  cmatrix_free(rho);

  return out;
}

cmatrix_t *vqe_noisy_prepare_density(int n_qubits, int n_layers,
                                     const double *theta, double gamma1,
                                     double gamma2, double gate_time) {
  if (n_qubits < 1 || n_layers < 1 || !theta || gamma1 < 0.0 || gamma2 < 0.0 ||
      gate_time < 0.0) {
    return NULL;
  }

  int dim = 1 << n_qubits;

  // \rho starts as |00...0><00...0|
  cmatrix_t *rho = cmatrix_alloc(dim, dim);
  if (!rho) {
    return NULL;
  }

  for (int i = 0; i < dim * dim; i++) {
    rho->data[i] = c_zero();
  }

  CMAT(rho, 0, 0) = c_real(1.0);

  // H=0 for the pure-decoherence steps between gates: only dissipators act, no
  // unitary drift beyond what the gates themselves already apply.
  cmatrix_t *Hzero = cmatrix_alloc(dim, dim);
  if (!Hzero) {
    cmatrix_free(rho);

    return NULL;
  }

  for (int i = 0; i < dim * dim; i++) {
    Hzero->data[i] = c_zero();
  }

  int n_ops = 2 * n_qubits;
  cmatrix_t **L = calloc((size_t)n_ops, sizeof *L);
  if (!L) {
    cmatrix_free(rho);
    cmatrix_free(Hzero);

    return NULL;
  }

  int alloc_ok = 1;
  for (int q = 0; q < n_qubits && alloc_ok; q++) {
    L[2 * q] = lindblad_amplitude_damping_op(n_qubits, q, gamma1);
    L[2 * q + 1] = lindblad_dephasing_op(n_qubits, q, gamma2);

    if (!L[2 * q] || !L[2 * q + 1]) {
      alloc_ok = 0;
    }
  }

  if (!alloc_ok) {
    for (int k = 0; k < n_ops; k++) {
      cmatrix_free(L[k]);
    }

    free(L);
    cmatrix_free(rho);
    cmatrix_free(Hzero);

    return NULL;
  }

  int ok = 1;
  for (int layer = 0; layer < n_layers && ok; layer++) {
    for (int q = 0; q < n_qubits && ok; q++) {
      complex_t g[4];
      ry_gate(theta[layer * n_qubits + q], g);

      cmatrix_t *U = embed_single_qubit_op(g, n_qubits, q);
      if (!U) {
        ok = 0;

        break;
      }

      rho = dm_conjugate_free(U, rho);
      cmatrix_free(U);
      if (!rho || lindblad_step_rk4(rho, Hzero, L, n_ops, gate_time) != 0) {
        ok = 0;

        break;
      }
    }

    for (int q = 0; q < n_qubits - 1 && ok; q++) {
      rho = dm_apply_cnot_free(rho, n_qubits, q, q + 1);
      if (!rho || lindblad_step_rk4(rho, Hzero, L, n_ops, gate_time) != 0) {
        ok = 0;

        break;
      }
    }
  }

  for (int k = 0; k < n_ops; k++) {
    cmatrix_free(L[k]);
  }
  free(L);
  cmatrix_free(Hzero);

  if (!ok) {
    cmatrix_free(rho);

    return NULL;
  }

  return rho;
}

double vqe_noisy_energy(int n_qubits, int n_layers, const double *theta,
                        const cmatrix_t *H, double gamma1, double gamma2,
                        double gate_time) {
  if (!H || n_qubits < 1 || H->nrows != H->ncols ||
      H->nrows != (1 << n_qubits)) {
    return 0.0;
  }

  cmatrix_t *rho = vqe_noisy_prepare_density(n_qubits, n_layers, theta, gamma1,
                                             gamma2, gate_time);
  if (!rho) {
    return 0.0;
  }

  // Tr(H \rho) = \sum_ij H[i][j] * \rho[j][i]
  complex_t trace = c_zero();
  int n = H->nrows;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      trace = c_add(trace, c_mul(CMAT(H, i, j), CMAT(rho, j, i)));
    }
  }

  cmatrix_free(rho);

  return trace.re;
}

/*
 * Solve the (d + 1)x(d + 1) normal-equations system for a degree-d
 * least-squares polynomial fit y ~= \sum_k p[k] * c^k through the given (c[i],
 * y[i]) points, via Gaussian elimination with partial pivoting. Writes p[0..d]
 * (p[0] is the c=0 intercept, to satisfy ZNE).
 *
 * Returns 0 on success, -1 on a singular system (shouldn't happen for distinct
 * c[i] with n_points >= d+1).
 */
static int polyfit_intercept(const double *c, const double *y, int n_points,
                             int degree, double *p_out) {
  int m = degree + 1;

  // Normal equations: (V^T V) p = V^T y, V[i][k] = c[i]^k
  double *ata = calloc((size_t)m * m, sizeof *ata);
  double *aty = calloc((size_t)m, sizeof *aty);
  double *cpow = malloc((size_t)m * sizeof *cpow);
  if (!ata || !aty || !cpow) {
    free(ata);
    free(aty);
    free(cpow);

    return -1;
  }

  for (int i = 0; i < n_points; i++) {
    cpow[0] = 1.0;
    for (int k = 1; k < m; k++) {
      cpow[k] = cpow[k - 1] * c[i];
    }

    for (int a = 0; a < m; a++) {
      aty[a] += cpow[a] * y[i];

      for (int b = 0; b < m; b++) {
        ata[a * m + b] += cpow[a] * cpow[b];
      }
    }
  }

  free(cpow);

  // Gaussian elimination with partial pivoting on [ata | aty]
  for (int col = 0; col < m; col++) {
    int piv = col;
    double best = fabs(ata[col * m + col]);

    for (int r = col + 1; r < m; r++) {
      double v = fabs(ata[r * m + col]);

      if (v > best) {
        best = v;
        piv = r;
      }
    }

    if (best < 1e-14) {
      free(ata);
      free(aty);

      return -1;
    }

    if (piv != col) {
      for (int k = 0; k < m; k++) {
        double tmp = ata[col * m + k];

        ata[col * m + k] = ata[piv * m + k];
        ata[piv * m + k] = tmp;
      }

      double tmp = aty[col];

      aty[col] = aty[piv];
      aty[piv] = tmp;
    }

    for (int r = col + 1; r < m; r++) {
      double factor = ata[r * m + col] / ata[col * m + col];

      for (int k = col; k < m; k++) {
        ata[r * m + k] -= factor * ata[col * m + k];
      }

      aty[r] -= factor * aty[col];
    }
  }

  for (int row = m - 1; row >= 0; row--) {
    double sum = aty[row];

    for (int k = row + 1; k < m; k++) {
      sum -= ata[row * m + k] * p_out[k];
    }

    p_out[row] = sum / ata[row * m + row];
  }

  free(ata);
  free(aty);

  return 0;
}

double vqe_noisy_zne_energy(int n_qubits, int n_layers, const double *theta,
                            const cmatrix_t *H, double gamma1, double gamma2,
                            double gate_time, int n_scales, double *raw_c1) {
  if (n_scales < 2) {
    return 0.0;
  }

  double *c = malloc((size_t)n_scales * sizeof *c);
  double *E = malloc((size_t)n_scales * sizeof *E);
  if (!c || !E) {
    free(c);
    free(E);

    return 0.0;
  }

  for (int k = 0; k < n_scales; k++) {
    c[k] = (double)(k + 1);
    E[k] = vqe_noisy_energy(n_qubits, n_layers, theta, H, gamma1 * c[k],
                            gamma2 * c[k], gate_time * c[k]);
  }

  if (raw_c1) {
    *raw_c1 = E[0];
  }

  int degree = n_scales - 1;
  double *p = malloc((size_t)(degree + 1) * sizeof *p);
  double result = E[0]; // fallback: no mitigation if the fit fails
  if (p && polyfit_intercept(c, E, n_scales, degree, p) == 0) {
    result = p[0];
  }

  free(p);
  free(c);
  free(E);

  return result;
}
