/*
Minimal multi-qubit state vector substrate
*/

#include "qubits.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "angular.h"
#include <math.h>
#include <stdlib.h>

const complex_t hadamard_gate[4] = {
    {M_SQRT1_2, 0.0}, {M_SQRT1_2, 0.0}, {M_SQRT1_2, 0.0}, {-M_SQRT1_2, 0.0}};

void rx_gate(double theta, complex_t out[4]) {
  double c = cos(theta / 2.0);
  double s = sin(theta / 2.0);

  out[0] = c_real(c);
  out[1] = c_new(0.0, -s);
  out[2] = c_new(0.0, -s);
  out[3] = c_real(c);
}

void ry_gate(double theta, complex_t out[4]) {
  double c = cos(theta / 2.0);
  double s = sin(theta / 2.0);

  out[0] = c_real(c);
  out[1] = c_real(-s);
  out[2] = c_real(s);
  out[3] = c_real(c);
}

void rz_gate(double theta, complex_t out[4]) {
  out[0] = c_new(cos(-theta / 2.0), sin(-theta / 2.0));
  out[1] = c_zero();
  out[2] = c_zero();
  out[3] = c_new(cos(theta / 2.0), sin(theta / 2.0));
}

cvector_t *qstate_alloc(int n_qubits) {
  if (n_qubits < 1) {
    return NULL;
  }

  long long dim = 1LL << n_qubits;
  cvector_t *psi = cvector_alloc((int)dim);
  if (!psi) {
    return NULL;
  }

  for (long long i = 0; i < dim; i++) {
    psi->data[i] = c_zero();
  }

  psi->data[0] = c_real(1.0);

  return psi;
}

double qstate_probability(const cvector_t *psi, int index) {
  if (!psi || index < 0 || index >= psi->n) {
    return 0.0;
  }

  return c_abs2(psi->data[index]);
}

void qstate_apply_gate1(cvector_t *psi, int n_qubits, int target,
                        const complex_t gate[4]) {
  if (!psi || target < 0 || target >= n_qubits) {
    return;
  }

  long long n = 1LL << n_qubits;
  int bitpos = n_qubits - 1 - target; // qubit 0 = leftmost/MSB
  long long mask = 1LL << bitpos;

  for (long long i = 0; i < n; i++) {
    if ((i & mask) == 0) {
      long long j = i | mask;
      complex_t a0 = psi->data[i], a1 = psi->data[j];
      psi->data[i] = c_add(c_mul(gate[0], a0), c_mul(gate[1], a1));
      psi->data[j] = c_add(c_mul(gate[2], a0), c_mul(gate[3], a1));
    }
  }
}

void qstate_apply_controlled_u(cvector_t *psi, int n_qubits, int control,
                               int target, const complex_t U[4]) {
  if (!psi || control < 0 || control >= n_qubits || target < 0 ||
      target >= n_qubits || control == target) {
    return;
  }

  long long n = 1LL << n_qubits;
  int cbit = n_qubits - 1 - control;
  int tbit = n_qubits - 1 - target;
  long long cmask = 1LL << cbit;
  long long tmask = 1LL << tbit;

  for (long long i = 0; i < n; i++) {
    // Only act where control=1 and, to touch each pair once, target=0
    if ((i & cmask) && (i & tmask) == 0) {
      long long j = i | tmask;
      complex_t a0 = psi->data[i], a1 = psi->data[j];
      psi->data[i] = c_add(c_mul(U[0], a0), c_mul(U[1], a1));
      psi->data[j] = c_add(c_mul(U[2], a0), c_mul(U[3], a1));
    }
  }
}

void qstate_apply_cnot(cvector_t *psi, int n_qubits, int control, int target) {
  qstate_apply_controlled_u(psi, n_qubits, control, target, sigma_x);
}

cmatrix_t *qstate_reduced_density_single(const cvector_t *psi, int n_qubits,
                                         int qubit) {
  if (!psi || qubit < 0 || qubit >= n_qubits) {
    return NULL;
  }

  cmatrix_t *rho = cmatrix_alloc(2, 2);
  if (!rho) {
    return NULL;
  }

  for (int a = 0; a < 2; a++) {
    for (int b = 0; b < 2; b++) {
      CMAT(rho, a, b) = c_zero();
    }
  }

  long long n = 1LL << n_qubits;
  int bitpos = n_qubits - 1 - qubit;
  long long mask = 1LL << bitpos;

  for (long long i = 0; i < n; i++) {
    if ((i & mask) != 0) {
      continue; // process each "rest" configuration once, at qubit=0
    }

    long long i0 = i;
    long long i1 = i | mask;
    complex_t p0 = psi->data[i0], p1 = psi->data[i1];

    CMAT(rho, 0, 0) = c_add(CMAT(rho, 0, 0), c_mul(c_conj(p0), p0));
    CMAT(rho, 0, 1) = c_add(CMAT(rho, 0, 1), c_mul(c_conj(p0), p1));
    CMAT(rho, 1, 0) = c_add(CMAT(rho, 1, 0), c_mul(c_conj(p1), p0));
    CMAT(rho, 1, 1) = c_add(CMAT(rho, 1, 1), c_mul(c_conj(p1), p1));
  }

  return rho;
}

double von_neumann_entropy_2x2(cmatrix_t *rho) {
  if (!rho) {
    return 0.0;
  }

  double a = CMAT(rho, 0, 0).re; // real: diagonal of a Hermitian matrix
  double d = CMAT(rho, 1, 1).re;
  complex_t b = CMAT(rho, 0, 1);
  double tr = a + d;
  double det = a * d - c_abs2(b); // real for Hermitian 2x2 (c = conj(b))

  double disc = tr * tr - 4.0 * det;
  if (disc < 0.0) {
    disc = 0.0; // numerical noise guard
  }

  double sq = sqrt(disc);
  double lambda1 = 0.5 * (tr + sq);
  double lambda2 = 0.5 * (tr - sq);

  double entropy = 0.0;
  double eigs[2] = {lambda1, lambda2};
  for (int k = 0; k < 2; k++) {
    if (eigs[k] > 1e-12) {
      entropy -= eigs[k] * log2(eigs[k]);
    }
  }

  return entropy;
}

int qstate_measure(cvector_t *psi, double u) {
  if (!psi || psi->n < 1) {
    return -1;
  }

  if (u < 0.0) {
    u = 0.0;
  }
  if (u >= 1.0) {
    u = 1.0 - 1e-15;
  }

  double cumulative = 0.0;
  int outcome = psi->n - 1; // fallback for floating-point round-off at u->1
  for (int i = 0; i < psi->n; i++) {
    cumulative += c_abs2(psi->data[i]);
    if (u < cumulative) {
      outcome = i;

      break;
    }
  }

  for (int i = 0; i < psi->n; i++) {
    psi->data[i] = c_zero();
  }

  psi->data[outcome] = c_real(1.0);

  return outcome;
}

int qstate_measure_qubit(cvector_t *psi, int n_qubits, int target, double u) {
  if (!psi || n_qubits < 1 || target < 0 || target >= n_qubits) {
    return -1;
  }

  long long n = 1LL << n_qubits;
  if (psi->n != (int)n) {
    return -1;
  }

  int bitpos = n_qubits - 1 - target;
  long long mask = 1LL << bitpos;

  double p0 = 0.0;
  for (long long i = 0; i < n; i++) {
    if ((i & mask) == 0) {
      p0 += c_abs2(psi->data[i]);
    }
  }

  if (u < 0.0) {
    u = 0.0;
  }
  if (u >= 1.0) {
    u = 1.0 - 1e-15;
  }

  int outcome = (u < p0) ? 0 : 1;
  double branch_prob = (outcome == 0) ? p0 : (1.0 - p0);

  double inv_norm = (branch_prob > 1e-300) ? 1.0 / sqrt(branch_prob) : 0.0;

  for (long long i = 0; i < n; i++) {
    int bit = (i & mask) != 0;
    if (bit == outcome) {
      psi->data[i] = c_scale(psi->data[i], inv_norm);
    } else {
      psi->data[i] = c_zero();
    }
  }

  return outcome;
}
