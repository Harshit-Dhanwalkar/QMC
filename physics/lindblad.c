/*
Lindblad (GKSL) equation: open quantum systems, density matrix evolution, and
single-qubit noise channels:

  d\rho/dt = -i[H, \rho] + \sum_k( L_k \rho L_k^\dagger - 1/2 {L_k^\dagger L_k,
\rho})
*/

#include "lindblad.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/ode/rk4.h"
#include "../core/vector.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  const cmatrix_t *H;
  cmatrix_t **L;
  int n_ops;
  int dim; // \rho is dim x dim
} lindblad_ctx_t;

cmatrix_t *density_from_pure_state(const cvector_t *psi) {
  if (!psi) {
    return NULL;
  }

  int n = psi->n;
  cmatrix_t *rho = cmatrix_alloc(n, n);
  if (!rho) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(rho, i, j) = c_mul(psi->data[i], c_conj(psi->data[j]));
    }
  }

  return rho;
}

double density_purity(const cmatrix_t *rho) {
  if (!rho) {
    return 0.0;
  }

  double sum = 0.0;

  for (int i = 0; i < rho->nrows * rho->ncols; i++) {
    sum += c_abs2(rho->data[i]);
  }

  return sum;
}

double density_von_neumann_entropy(const cmatrix_t *rho) {
  if (!rho) {
    return 0.0;
  }

  eigen_t *e = cmatrix_eigh(rho);
  if (!e) {
    return 0.0;
  }

  double entropy = 0.0;
  for (int i = 0; i < e->n; i++) {
    double lambda = e->eigenvalues[i];
    if (lambda > 1e-12) {
      entropy -= lambda * log2(lambda);
    }
  }

  eigen_free(e);
  return entropy;
}

static void lindblad_ode_rhs(double t, const cvector_t *y, cvector_t *dydt,
                             void *params) {
  lindblad_ctx_t *ctx = (lindblad_ctx_t *)params;
  int n = ctx->dim;

  cmatrix_t rho_view = {.data = y->data, .nrows = n, .ncols = n};
  cmatrix_t *result = lindblad_rhs(ctx->H, &rho_view, ctx->L, ctx->n_ops);

  if (!result) {
    for (int i = 0; i < n * n; i++) {
      dydt->data[i] = c_zero();
    }
    return;
  }

  memcpy(dydt->data, result->data, (size_t)n * n * sizeof(complex_t));
  cmatrix_free(result);
}

cmatrix_t *lindblad_rhs(const cmatrix_t *H, const cmatrix_t *rho, cmatrix_t **L,
                        int n_ops) {
  if (!H || !rho) {
    return NULL;
  }

  int n = rho->nrows;
  if (rho->ncols != n || H->nrows != n || H->ncols != n) {
    return NULL;
  }

  cmatrix_t *Hrho = cmatrix_multiply(H, rho);
  cmatrix_t *rhoH = cmatrix_multiply(rho, H);
  if (!Hrho || !rhoH) {
    cmatrix_free(Hrho);
    cmatrix_free(rhoH);

    return NULL;
  }

  cmatrix_t *result = cmatrix_alloc(n, n);
  if (!result) {
    cmatrix_free(Hrho);
    cmatrix_free(rhoH);

    return NULL;
  }

  // -i[H, \rho] = -i * (H * \rho - \rho * H)
  for (int i = 0; i < n * n; i++) {
    complex_t comm = c_sub(Hrho->data[i], rhoH->data[i]);

    result->data[i] = c_mul(c_imag(-1.0), comm);
  }

  cmatrix_free(Hrho);
  cmatrix_free(rhoH);

  for (int k = 0; k < n_ops; k++) {
    const cmatrix_t *Lk = L ? L[k] : NULL;
    if (!Lk || Lk->nrows != n || Lk->ncols != n) {
      cmatrix_free(result);

      return NULL;
    }

    cmatrix_t *Ldag = cmatrix_adjoint(Lk);
    cmatrix_t *LdagL = Ldag ? cmatrix_multiply(Ldag, Lk) : NULL;
    cmatrix_t *Lrho = cmatrix_multiply(Lk, rho);
    cmatrix_t *LrhoLdag = (Lrho && Ldag) ? cmatrix_multiply(Lrho, Ldag) : NULL;
    cmatrix_t *LdagLrho = LdagL ? cmatrix_multiply(LdagL, rho) : NULL;
    cmatrix_t *rhoLdagL = LdagL ? cmatrix_multiply(rho, LdagL) : NULL;

    if (!Ldag || !LdagL || !Lrho || !LrhoLdag || !LdagLrho || !rhoLdagL) {
      cmatrix_free(Ldag);
      cmatrix_free(LdagL);
      cmatrix_free(Lrho);
      cmatrix_free(LrhoLdag);
      cmatrix_free(LdagLrho);
      cmatrix_free(rhoLdagL);
      cmatrix_free(result);

      return NULL;
    }

    for (int i = 0; i < n * n; i++) {
      complex_t anticomm = c_add(LdagLrho->data[i], rhoLdagL->data[i]);
      complex_t dissipator = c_sub(LrhoLdag->data[i], c_scale(anticomm, 0.5));

      result->data[i] = c_add(result->data[i], dissipator);
    }

    cmatrix_free(Ldag);
    cmatrix_free(LdagL);
    cmatrix_free(Lrho);
    cmatrix_free(LrhoLdag);
    cmatrix_free(LdagLrho);
    cmatrix_free(rhoLdagL);
  }

  return result;
}

int lindblad_step_rk4(const cmatrix_t *rho, const cmatrix_t *H, cmatrix_t **L,
                      int n_ops, double dt) {
  if (!rho || !H) {
    return -1;
  }

  int n = rho->nrows;
  if (rho->ncols != n || H->nrows != n || H->ncols != n) {
    return -1;
  }

  cvector_t rho_flat = {.data = rho->data, .n = n * n};
  lindblad_ctx_t ctx = {H, L, n_ops, n};

  // cmatrix_t *tmp = cmatrix_alloc(n, n);
  // if (!tmp) {
  //   return -1;
  // }
  //
  // cmatrix_t *k1 = lindblad_rhs(H, rho, L, n_ops);
  // if (!k1) {
  //   cmatrix_free(tmp);
  //   return -1;
  // }
  //
  // for (int i = 0; i < n * n; i++) {
  //   tmp->data[i] = c_add(rho->data[i], c_scale(k1->data[i], dt * 0.5));
  // }
  //
  // cmatrix_t *k2 = lindblad_rhs(H, tmp, L, n_ops);
  // if (!k2) {
  //   cmatrix_free(k1);
  //   cmatrix_free(tmp);
  //
  //   return -1;
  // }
  //
  // for (int i = 0; i < n * n; i++) {
  //   tmp->data[i] = c_add(rho->data[i], c_scale(k2->data[i], dt * 0.5));
  // }
  //
  // cmatrix_t *k3 = lindblad_rhs(H, tmp, L, n_ops);
  // if (!k3) {
  //   cmatrix_free(k1);
  //   cmatrix_free(k2);
  //   cmatrix_free(tmp);
  //
  //   return -1;
  // }
  //
  // for (int i = 0; i < n * n; i++) {
  //   tmp->data[i] = c_add(rho->data[i], c_scale(k3->data[i], dt));
  // }
  //
  // cmatrix_t *k4 = lindblad_rhs(H, tmp, L, n_ops);
  // if (!k4) {
  //   cmatrix_free(k1);
  //   cmatrix_free(k2);
  //   cmatrix_free(k3);
  //   cmatrix_free(tmp);
  //
  //   return -1;
  // }
  //
  // for (int i = 0; i < n * n; i++) {
  //   complex_t sum =
  //       c_add(k1->data[i], c_scale(c_add(k2->data[i], k3->data[i]), 2.0));
  //   sum = c_add(sum, k4->data[i]);
  //   rho->data[i] = c_add(rho->data[i], c_scale(sum, dt / 6.0));
  // }
  //
  // cmatrix_free(k1);
  // cmatrix_free(k2);
  // cmatrix_free(k3);
  // cmatrix_free(k4);
  // cmatrix_free(tmp);

  return rk4_step(0.0, dt, &rho_flat, lindblad_ode_rhs, &ctx);
}

int lindblad_evolve(const cmatrix_t *rho, const cmatrix_t *H, cmatrix_t **L,
                    int n_ops, double dt, int steps) {
  if (!rho || !H || steps < 0) {
    return -1;
  }

  for (int s = 0; s < steps; s++) {
    if (lindblad_step_rk4(rho, H, L, n_ops, dt) != 0)
      return -1;
  }

  return 0;
}

cmatrix_t *embed_single_qubit_op(const complex_t op[4], int n_qubits,
                                 int target) {
  if (!op || n_qubits < 1 || target < 0 || target >= n_qubits) {
    return NULL;
  }

  int dim = 1 << n_qubits;
  cmatrix_t *M = cmatrix_alloc(dim, dim);
  if (!M) {
    return NULL;
  }

  int bitpos = n_qubits - 1 - target; // qubit 0 = leftmost/MSB
  int mask = 1 << bitpos;

  for (int r = 0; r < dim; r++) {
    for (int c = 0; c < dim; c++) {
      int r_rest = r & ~mask;
      int c_rest = c & ~mask;

      if (r_rest != c_rest) {
        CMAT(M, r, c) = c_zero();

        continue;
      }

      int a = (r & mask) ? 1 : 0; // target-qubit bit of row
      int b = (c & mask) ? 1 : 0; // target-qubit bit of column
      CMAT(M, r, c) = op[a * 2 + b];
    }
  }

  return M;
}

cmatrix_t *lindblad_amplitude_damping_op(int n_qubits, int target,
                                         double gamma) {
  // sigma_minus = |0><1| : lowers |1> -> |0>
  const complex_t sigma_minus[4] = {
      {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  cmatrix_t *M = embed_single_qubit_op(sigma_minus, n_qubits, target);

  if (!M) {
    return NULL;
  }

  cmatrix_scale(M, c_real(sqrt(gamma)));

  return M;
}

cmatrix_t *lindblad_dephasing_op(int n_qubits, int target, double gamma) {
  const complex_t sigma_z[4] = {
      {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {-1.0, 0.0}};
  cmatrix_t *M = embed_single_qubit_op(sigma_z, n_qubits, target);

  if (!M) {
    return NULL;
  }

  cmatrix_scale(M, c_real(sqrt(gamma / 2.0)));

  return M;
}

cmatrix_t *lindblad_bitflip_op(int n_qubits, int target, double gamma) {
  const complex_t sigma_x[4] = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}};
  cmatrix_t *M = embed_single_qubit_op(sigma_x, n_qubits, target);

  if (!M) {
    return NULL;
  }

  cmatrix_scale(M, c_real(sqrt(gamma)));

  return M;
}

int density_measure_computational_basis(cmatrix_t *rho, double u) {
  if (!rho || rho->nrows != rho->ncols) {
    return -1;
  }

  int n = rho->nrows;

  if (u < 0.0) {
    u = 0.0;
  }
  if (u >= 1.0) {
    u = 1.0 - 1e-15;
  }

  double cumulative = 0.0;
  int outcome = n - 1; // fallback for floating-point round-off at u -> 1
  for (int i = 0; i < n; i++) {
    double p = CMAT(rho, i, i).re;
    cumulative += p;
    if (u < cumulative) {
      outcome = i;

      break;
    }
  }

  for (int i = 0; i < n * n; i++) {
    rho->data[i] = c_zero();
  }

  CMAT(rho, outcome, outcome) = c_real(1.0);

  return outcome;
}
