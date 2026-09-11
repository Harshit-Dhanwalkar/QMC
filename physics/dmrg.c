/*
 * Infinite-system DMRG for the open-boundary spin-1/2 XXZ chain
 */

#include "dmrg.h"
#include "../core/complex.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include <math.h>
#include <stdlib.h>

/* Bare 2x2 site operators basis order [up, down], S^z_up=+1/2, S^z_down=-1/2,
 * S^+ raises down->up. */
static cmatrix_t *site_sz(void) {
  cmatrix_t *m = cmatrix_alloc(2, 2);
  if (!m) {
    return NULL;
  }

  CMAT(m, 0, 0) = c_real(0.5);
  CMAT(m, 0, 1) = c_zero();
  CMAT(m, 1, 0) = c_zero();
  CMAT(m, 1, 1) = c_real(-0.5);

  return m;
}

static cmatrix_t *site_sp(void) {
  cmatrix_t *m = cmatrix_alloc(2, 2);
  if (!m) {
    return NULL;
  }

  CMAT(m, 0, 0) = c_zero();
  CMAT(m, 0, 1) = c_real(1.0); // |up><down|
  CMAT(m, 1, 0) = c_zero();
  CMAT(m, 1, 1) = c_zero();

  return m;
}

static cmatrix_t *identity(int n) {
  cmatrix_t *m = cmatrix_alloc(n, n);
  if (!m) {
    return NULL;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      CMAT(m, i, j) = (i == j) ? c_real(1.0) : c_zero();
    }
  }

  return m;
}

dmrg_block_t *dmrg_block_init(void) {
  dmrg_block_t *b = malloc(sizeof(dmrg_block_t));
  if (!b) {
    return NULL;
  }

  b->dim = 2;
  b->H = cmatrix_alloc(2, 2);
  b->Sz_end = site_sz();
  b->Sp_end = site_sp();
  if (!b->H || !b->Sz_end || !b->Sp_end) {
    dmrg_block_free(b);

    return NULL;
  }

  for (int i = 0; i < 4; i++) {
    b->H->data[i] = c_zero();
  }

  return b;
}

void dmrg_block_free(dmrg_block_t *b) {
  if (!b) {
    return;
  }

  cmatrix_free(b->H);
  cmatrix_free(b->Sz_end);
  cmatrix_free(b->Sp_end);
  free(b);
}

/*
 * out = a + scal e *b (element-wise), in a freshly allocated matrix
 *
 * Returns NULL on allocation failure.
 */
static cmatrix_t *add_scaled(const cmatrix_t *a, const cmatrix_t *b,
                             double scale) {
  cmatrix_t *bs = cmatrix_copy(b);
  if (!bs) {
    return NULL;
  }

  cmatrix_scale(bs, c_real(scale));
  cmatrix_t *out = cmatrix_add(a, bs);
  cmatrix_free(bs);

  return out;
}

/*
 * Enlarge block `b` by one bare site: new dim = b->dim * 2
 * H_new  = H_b (x) I2 + Jz*(Sz_b (x) Sz) + (Jxy/2)*(Sp_b (x) Sm + Sm_b (x) Sp)
 * Sz_new = I_dim (x) Sz   (acts on the newly appended site)
 * Sp_new = I_dim (x) Sp
 *
 * Returns NULL on allocation failure; frees no inputs
 */
static dmrg_block_t *enlarge_block(const dmrg_block_t *b, double Jz, double Jxy,
                                   const cmatrix_t *sz, const cmatrix_t *sp,
                                   const cmatrix_t *sm) {
  dmrg_block_t *nb = malloc(sizeof(dmrg_block_t));
  if (!nb) {
    return NULL;
  }

  nb->dim = b->dim * 2;
  nb->H = NULL;
  nb->Sz_end = NULL;
  nb->Sp_end = NULL;

  cmatrix_t *I2 = identity(2);
  cmatrix_t *I_dim = identity(b->dim);
  cmatrix_t *sm_b = cmatrix_adjoint(b->Sp_end); // S-_b = (S+_b)^\dagger
  if (!I2 || !I_dim || !sm_b) {
    goto fail;
  }

  cmatrix_t *term1 = cmatrix_kron(b->H, I2);
  cmatrix_t *term2 = cmatrix_kron(b->Sz_end, sz);
  cmatrix_t *term3 = cmatrix_kron(b->Sp_end, sm);
  cmatrix_t *term4 = cmatrix_kron(sm_b, sp);
  if (!term1 || !term2 || !term3 || !term4) {
    cmatrix_free(term1);
    cmatrix_free(term2);
    cmatrix_free(term3);
    cmatrix_free(term4);
    goto fail;
  }

  cmatrix_t *h1 = add_scaled(term1, term2, Jz);
  cmatrix_t *h2 = add_scaled(term3, term4, 1.0);
  cmatrix_free(term1);
  cmatrix_free(term2);
  cmatrix_free(term3);
  cmatrix_free(term4);
  if (!h1 || !h2) {
    cmatrix_free(h1);
    cmatrix_free(h2);
    goto fail;
  }

  nb->H = add_scaled(h1, h2, Jxy / 2.0);
  cmatrix_free(h1);
  cmatrix_free(h2);

  nb->Sz_end = cmatrix_kron(I_dim, sz);
  nb->Sp_end = cmatrix_kron(I_dim, sp);

  cmatrix_free(I2);
  cmatrix_free(I_dim);
  cmatrix_free(sm_b);

  if (!nb->H || !nb->Sz_end || !nb->Sp_end) {
    dmrg_block_free(nb);

    return NULL;
  }

  return nb;

fail:
  cmatrix_free(I2);
  cmatrix_free(I_dim);
  cmatrix_free(sm_b);
  free(nb);

  return NULL;
}

/*
 * Matrix-free application of the superblock Hamiltonian to V (D x D, system
 * index x mirrored-environment index):
 *  H_super = H_eb(x)I + I(x)H_eb + Jz * (Sz_eb(x)Sz_eb)
 *                     + (Jxy / 2) * (Sp_eb(x)Sm_eb + Sm_eb(x)Sp_eb)
 *
 * Using the identity (A(x)B) applied to V (reshaped as a matrix) equals A @ V @
 * B^T, and that H_eb, Sz_eb are real-symmetric (B^T=B) while Sm_eb = Sp_eb^T
 * (so both cross terms reduce to Sp_eb@V@Sp_eb and Sm_eb@V@Sm_eb):
 *   H_super(V) = H_eb@V + V@H_eb + Jz*(Sz_eb@V@Sz_eb)
 *                       + (Jxy / 2)*(Sp_eb@V@Sp_eb + Sm_eb@V@Sm_eb)
 *
 * Returns NULL on allocation failure
 */
static cmatrix_t *apply_superblock(const dmrg_block_t *eb, double Jz,
                                   double Jxy, const cmatrix_t *V) {
  cmatrix_t *sm_eb = cmatrix_adjoint(eb->Sp_end);
  if (!sm_eb) {
    return NULL;
  }

  cmatrix_t *HV = cmatrix_multiply(eb->H, V);
  cmatrix_t *VH = cmatrix_multiply(V, eb->H);
  cmatrix_t *SzV = cmatrix_multiply(eb->Sz_end, V);
  cmatrix_t *SzVSz = SzV ? cmatrix_multiply(SzV, eb->Sz_end) : NULL;
  cmatrix_t *SpV = cmatrix_multiply(eb->Sp_end, V);
  cmatrix_t *SpVSp = SpV ? cmatrix_multiply(SpV, eb->Sp_end) : NULL;
  cmatrix_t *SmV = cmatrix_multiply(sm_eb, V);
  cmatrix_t *SmVSm = SmV ? cmatrix_multiply(SmV, sm_eb) : NULL;

  cmatrix_free(sm_eb);
  cmatrix_free(SzV);
  cmatrix_free(SpV);
  cmatrix_free(SmV);

  cmatrix_t *out = NULL;
  if (HV && VH && SzVSz && SpVSp && SmVSm) {
    cmatrix_t *s1 = cmatrix_add(HV, VH);
    cmatrix_t *s2 = s1 ? add_scaled(s1, SzVSz, Jz) : NULL;
    cmatrix_t *hop = cmatrix_add(SpVSp, SmVSm);

    if (s2 && hop) {
      out = add_scaled(s2, hop, Jxy / 2.0);
    }

    cmatrix_free(s1);
    cmatrix_free(s2);
    cmatrix_free(hop);
  }

  cmatrix_free(HV);
  cmatrix_free(VH);
  cmatrix_free(SzVSz);
  cmatrix_free(SpVSp);
  cmatrix_free(SmVSm);

  return out;
}

/*
 * Real part of Frobenius inner product <a|b> = \sum conj(a_i) b_i, treating a D
 * x D cmatrix_t as a flat vector of length D*D. Operators are all real-valued,
 * Lanczos recursion needs
 */
static double frob_dot_real(const cmatrix_t *a, const cmatrix_t *b) {
  int n = a->nrows * a->ncols;
  complex_t sum = c_zero();

  for (int i = 0; i < n; i++) {
    sum = c_add(sum, c_mul(c_conj(a->data[i]), b->data[i]));
  }

  return sum.re;
}

static double frob_norm(const cmatrix_t *a) {
  return sqrt(frob_dot_real(a, a));
}

/* y += \alpha * x, element-wise (\alpha real)*/
static void axpy_inplace(cmatrix_t *y, double alpha, const cmatrix_t *x) {
  int n = y->nrows * y->ncols;

  for (int i = 0; i < n; i++) {
    y->data[i] = c_add(y->data[i], c_scale(x->data[i], alpha));
  }
}

/*
 * Ground state of superblock via Lanczos with full reorthogonalization.
 *
 * Returns normalized ground-state D x D matrix and sets *E0_out, or NULL on
 * allocation/breakdown failure
 */
static cmatrix_t *lanczos_ground_state(const dmrg_block_t *eb, double Jz,
                                       double Jxy, double *E0_out) {
  int D = eb->dim;
  int dim = D * D;
  int max_iter = (dim < 80) ? dim : 80; // Krylov subspace <= problem dim
  const double tol = 1e-12;

  cmatrix_t **basis = calloc((size_t)max_iter + 1, sizeof(cmatrix_t *));
  double *alpha = malloc(sizeof(double) * (size_t)max_iter);
  double *beta = malloc(sizeof(double) * (size_t)max_iter);
  if (!basis || !alpha || !beta) {
    free(basis);
    free(alpha);
    free(beta);

    return NULL;
  }

  cmatrix_t *v0 = cmatrix_alloc(D, D);
  if (!v0) {
    free(basis);
    free(alpha);
    free(beta);

    return NULL;
  }

  for (int i = 0; i < dim; i++) {
    v0->data[i] = c_real(1.0);
  }

  double nrm0 = frob_norm(v0);
  for (int i = 0; i < dim; i++) {
    v0->data[i] = c_scale(v0->data[i], 1.0 / nrm0);
  }

  basis[0] = v0;

  int M = 0; // number of basis vectors / alpha coefficients generated
  for (int j = 0; j < max_iter; j++) {
    cmatrix_t *w = apply_superblock(eb, Jz, Jxy, basis[j]);
    if (!w) {
      goto lanczos_fail;
    }

    alpha[j] = frob_dot_real(basis[j], w);
    axpy_inplace(w, -alpha[j], basis[j]);
    if (j > 0) {
      axpy_inplace(w, -beta[j - 1], basis[j - 1]);
    }

    // full reorthogonalization against every prior basis vector
    for (int k = 0; k <= j; k++) {
      double proj = frob_dot_real(basis[k], w);

      axpy_inplace(w, -proj, basis[k]);
    }

    double bj = frob_norm(w);
    M = j + 1;
    if (bj < tol || j == max_iter - 1) {
      cmatrix_free(w);

      break;
    }

    beta[j] = bj;
    for (int i = 0; i < dim; i++) {
      w->data[i] = c_scale(w->data[i], 1.0 / bj);
    }

    basis[j + 1] = w;
  }

  eigen_t *teig = tridiag_eigh(alpha, beta, M); // ascending eigenvalues
  if (!teig) {
    goto lanczos_fail;
  }

  cmatrix_t *gs = cmatrix_alloc(D, D);
  if (!gs) {
    eigen_free(teig);

    goto lanczos_fail;
  }

  for (int i = 0; i < dim; i++) {
    gs->data[i] = c_zero();
  }

  for (int k = 0; k < M; k++) {
    double coeff = CMAT(teig->eigenvectors, k, 0).re;

    axpy_inplace(gs, coeff, basis[k]);
  }

  double gs_norm = frob_norm(gs);
  for (int i = 0; i < dim; i++) {
    gs->data[i] = c_scale(gs->data[i], 1.0 / gs_norm);
  }

  *E0_out = teig->eigenvalues[0];

  eigen_free(teig);

  for (int k = 0; k <= max_iter; k++) {
    cmatrix_free(basis[k]);
  }
  free(basis);
  free(alpha);
  free(beta);

  return gs;

lanczos_fail:
  for (int k = 0; k <= max_iter; k++) {
    cmatrix_free(basis[k]);
  }
  free(basis);
  free(alpha);
  free(beta);

  return NULL;
}

dmrg_result_t *dmrg_run(int N_target, double Jz, double Jxy, int m_max) {
  if (N_target < 2 || m_max < 1 || !isfinite(Jz) || !isfinite(Jxy)) {
    return NULL;
  }

  cmatrix_t *sz = site_sz();
  cmatrix_t *sp = site_sp();
  cmatrix_t *sm = sp ? cmatrix_adjoint(sp) : NULL;
  dmrg_block_t *block = dmrg_block_init();
  if (!sz || !sp || !sm || !block) {
    cmatrix_free(sz);
    cmatrix_free(sp);
    cmatrix_free(sm);
    dmrg_block_free(block);

    return NULL;
  }

  dmrg_result_t *result = malloc(sizeof(dmrg_result_t));
  if (!result) {
    cmatrix_free(sz);
    cmatrix_free(sp);
    cmatrix_free(sm);
    dmrg_block_free(block);

    return NULL;
  }

  result->truncation_dim = block->dim;
  result->truncation_error = 0.0;

  int L = 1; // sites currently in `block`

  for (;;) {
    dmrg_block_t *eb = enlarge_block(block, Jz, Jxy, sz, sp, sm);
    dmrg_block_free(block);
    block = NULL;
    if (!eb) {
      goto fail;
    }

    int L_new = L + 1;
    int total_N = 2 * L_new;
    int D = eb->dim;

    double E0;
    cmatrix_t *Psi = lanczos_ground_state(eb, Jz, Jxy, &E0);
    if (!Psi) {
      dmrg_block_free(eb);

      goto fail;
    }

    if (total_N >= N_target) {
      result->energy = E0;
      result->N_reached = total_N;
      result->energy_per_site = E0 / (double)total_N;

      cmatrix_free(Psi);
      dmrg_block_free(eb);
      cmatrix_free(sz);
      cmatrix_free(sp);
      cmatrix_free(sm);

      return result;
    }

    /* Reduced density matrix of the system (enlarged block) side:
     * Psi[a][e] (a=system index, e=mirrored-environment index) is ground state
     * reshaped as a D x D matrix; \rho = \Psi * \Psi^\dagger. */
    cmatrix_t *Psi_dag = cmatrix_adjoint(Psi);
    cmatrix_t *rho = Psi_dag ? cmatrix_multiply(Psi, Psi_dag) : NULL;

    cmatrix_free(Psi);
    cmatrix_free(Psi_dag);
    if (!rho) {
      dmrg_block_free(eb);

      goto fail;
    }

    eigen_t *rho_eig = cmatrix_eigh(rho); // ascending eigenvalues
    cmatrix_free(rho);
    if (!rho_eig) {
      dmrg_block_free(eb);

      goto fail;
    }

    int m = (D < m_max) ? D : m_max;
    cmatrix_t *O = cmatrix_alloc(D, m);
    if (!O) {
      eigen_free(rho_eig);
      dmrg_block_free(eb);

      goto fail;
    }

    double kept = 0.0;
    for (int col = 0; col < m; col++) {
      int src_col = D - 1 - col; // largest eigenvalues are at the end
      kept += rho_eig->eigenvalues[src_col];
      for (int row = 0; row < D; row++) {
        CMAT(O, row, col) = CMAT(rho_eig->eigenvectors, row, src_col);
      }
    }

    eigen_free(rho_eig);

    result->truncation_dim = m;
    result->truncation_error = 1.0 - kept;

    cmatrix_t *O_dag = cmatrix_adjoint(O);
    cmatrix_t *H_trunc = NULL, *Sz_trunc = NULL, *Sp_trunc = NULL;
    if (O_dag) {
      cmatrix_t *tmp;

      tmp = cmatrix_multiply(eb->H, O);
      H_trunc = tmp ? cmatrix_multiply(O_dag, tmp) : NULL;

      cmatrix_free(tmp);

      tmp = cmatrix_multiply(eb->Sz_end, O);
      Sz_trunc = tmp ? cmatrix_multiply(O_dag, tmp) : NULL;

      cmatrix_free(tmp);

      tmp = cmatrix_multiply(eb->Sp_end, O);
      Sp_trunc = tmp ? cmatrix_multiply(O_dag, tmp) : NULL;

      cmatrix_free(tmp);
    }

    cmatrix_free(O);
    cmatrix_free(O_dag);
    dmrg_block_free(eb);

    if (!H_trunc || !Sz_trunc || !Sp_trunc) {
      cmatrix_free(H_trunc);
      cmatrix_free(Sz_trunc);
      cmatrix_free(Sp_trunc);

      goto fail;
    }

    block = malloc(sizeof(dmrg_block_t));
    if (!block) {
      cmatrix_free(H_trunc);
      cmatrix_free(Sz_trunc);
      cmatrix_free(Sp_trunc);

      goto fail;
    }

    block->dim = m;
    block->H = H_trunc;
    block->Sz_end = Sz_trunc;
    block->Sp_end = Sp_trunc;
    L = L_new;
  }

fail:
  cmatrix_free(sz);
  cmatrix_free(sp);
  cmatrix_free(sm);
  dmrg_block_free(block);
  free(result);

  return NULL;
}
