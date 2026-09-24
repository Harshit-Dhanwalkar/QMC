/*
 * Finite-system DMRG sweeps for the open-boundary spin-1/2 XXZ chain
 *
 * This file is deliberately self-contained (it does not reach into dmrg.c's
 * static helpers) so validated infinite-system algorithm in dmrg.c is
 * untouched. The only things shared with dmrg.h are dmrg_block_t layout and
 * dmrg_block_init()/dmrg_block_free()
 */

#include "finite_dmrg.h"
#include "../core/complex.h"
#include "../core/linalg/tridiag_eigh.h"
#include "../core/matrix.h"
#include "physics/dmrg.h"
#include <math.h>
#include <stdlib.h>

/* Bare 2x2 site operators, basis order [up, down], S^z_up=+1/2,
 * S^z_down=-1/2, S^+ raises down->up */
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

// out = a + scale*b (element-wise), in a freshly allocated matrix
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
 * Enlarge block `b` by one bare site: new dim = b->dim * 2. Identical
 * construction to dmrg.c's enlarge_block. Returns NULL on allocation
 * failure; frees none of its inputs
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
  cmatrix_t *sm_b = cmatrix_adjoint(b->Sp_end);
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
 * Matrix-free application of a superblock Hamiltonian built from two distinct
 * blocks (system `sys`, environment `env`, dims D1 and D2) to V (D1 x D2,
 * system index x environment index):
 *
 *   H_super = H_sys(x)I_env + I_sys(x)H_env + Jz * (Sz_sys(x)Sz_env)
 *                  + (Jxy/2) * (Sp_sys(x)Sm_env + Sm_sys(x)Sp_env)
 *
 * Using (A(x)B) applied to V (reshaped as a matrix) = A @ V @ B^T, and that
 * H_env, Sz_env are real-symmetric (B^T=B) while Sm_env = Sp_env^T:
 *
 *   H_super(V) = H_sys@V + V@H_env + Jz * (Sz_sys@V@Sz_env)
 *                        + (Jxy/2) * (Sp_sys@V@Sp_env + Sm_sys@V@Sm_env)
 *
 * Returns NULL on allocation failure
 */
static cmatrix_t *apply_superblock_general(const dmrg_block_t *sys,
                                           const dmrg_block_t *env, double Jz,
                                           double Jxy, const cmatrix_t *V) {
  cmatrix_t *sm_sys = cmatrix_adjoint(sys->Sp_end);
  cmatrix_t *sm_env = cmatrix_adjoint(env->Sp_end);
  if (!sm_sys || !sm_env) {
    cmatrix_free(sm_sys);
    cmatrix_free(sm_env);

    return NULL;
  }

  cmatrix_t *HV = cmatrix_multiply(sys->H, V);
  cmatrix_t *VH = cmatrix_multiply(V, env->H);
  cmatrix_t *SzV = cmatrix_multiply(sys->Sz_end, V);
  cmatrix_t *SzVSz = SzV ? cmatrix_multiply(SzV, env->Sz_end) : NULL;
  cmatrix_t *SpV = cmatrix_multiply(sys->Sp_end, V);
  cmatrix_t *SpVSp = SpV ? cmatrix_multiply(SpV, env->Sp_end) : NULL;
  cmatrix_t *SmV = cmatrix_multiply(sm_sys, V);
  cmatrix_t *SmVSm = SmV ? cmatrix_multiply(SmV, sm_env) : NULL;

  cmatrix_free(sm_sys);
  cmatrix_free(sm_env);
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

// Real part of Frobenius inner product, treating a matrix as a flat vector
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

// y += \alpha * x, element-wise (\alpha real)
static void axpy_inplace(cmatrix_t *y, double alpha, const cmatrix_t *x) {
  int n = y->nrows * y->ncols;

  for (int i = 0; i < n; i++) {
    y->data[i] = c_add(y->data[i], c_scale(x->data[i], alpha));
  }
}

/*
 * Ground state of (sys, env) superblock via Lanczos with full
 * reorthogonalization
 * NOTE: Generalized from dmrg.c's lanczos_ground_state to allow sys->dim !=
 * env->dim
 *
 * Returns a normalized D1 x D2 ground-state matrix and sets *E0_out, or NULL on
 * allocation/breakdown failure
 */
static cmatrix_t *lanczos_ground_state_general(const dmrg_block_t *sys,
                                               const dmrg_block_t *env,
                                               double Jz, double Jxy,
                                               double *E0_out) {
  int D1 = sys->dim;
  int D2 = env->dim;
  int dim = D1 * D2;
  int max_iter = (dim < 80) ? dim : 80;
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

  cmatrix_t *v0 = cmatrix_alloc(D1, D2);
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

  int M = 0;
  for (int j = 0; j < max_iter; j++) {
    cmatrix_t *w = apply_superblock_general(sys, env, Jz, Jxy, basis[j]);
    if (!w) {
      goto lanczos_fail;
    }

    alpha[j] = frob_dot_real(basis[j], w);
    axpy_inplace(w, -alpha[j], basis[j]);
    if (j > 0) {
      axpy_inplace(w, -beta[j - 1], basis[j - 1]);
    }

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

  eigen_t *teig = tridiag_eigh(alpha, beta, M);
  if (!teig) {
    goto lanczos_fail;
  }

  cmatrix_t *gs = cmatrix_alloc(D1, D2);
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

/*
 * Truncate D x D `raw` block (just-enlarged system side) down to at most m_max
 * states, using the reduced density matrix of ground state \Psi (D x D2, D2 =
 * environment dim) on raw/system (row) side:
 *   \rho = \Psi @ \Psi^\dagger
 *
 * Returns a freshly allocated truncated dmrg_block_t and sets *trunc_err_out,
 * or NULL on allocation failure
 */
static dmrg_block_t *truncate_block(const dmrg_block_t *raw,
                                    const cmatrix_t *Psi, int m_max,
                                    double *trunc_err_out) {
  cmatrix_t *Psi_dag = cmatrix_adjoint(Psi);
  cmatrix_t *rho = Psi_dag ? cmatrix_multiply(Psi, Psi_dag) : NULL;

  cmatrix_free(Psi_dag);
  if (!rho) {
    return NULL;
  }

  eigen_t *rho_eig = cmatrix_eigh(rho); // ascending eigenvalues
  cmatrix_free(rho);
  if (!rho_eig) {
    return NULL;
  }

  int D = raw->dim;
  int m = (D < m_max) ? D : m_max;
  cmatrix_t *O = cmatrix_alloc(D, m);
  if (!O) {
    eigen_free(rho_eig);

    return NULL;
  }

  double kept = 0.0;
  for (int col = 0; col < m; col++) {
    int src_col = D - 1 - col; // largest eigenvalues are at end

    kept += rho_eig->eigenvalues[src_col];
    for (int row = 0; row < D; row++) {
      CMAT(O, row, col) = CMAT(rho_eig->eigenvectors, row, src_col);
    }
  }

  eigen_free(rho_eig);

  *trunc_err_out = 1.0 - kept;

  cmatrix_t *O_dag = cmatrix_adjoint(O);
  cmatrix_t *H_trunc = NULL, *Sz_trunc = NULL, *Sp_trunc = NULL;
  if (O_dag) {
    cmatrix_t *tmp;

    tmp = cmatrix_multiply(raw->H, O);
    H_trunc = tmp ? cmatrix_multiply(O_dag, tmp) : NULL;

    cmatrix_free(tmp);

    tmp = cmatrix_multiply(raw->Sz_end, O);
    Sz_trunc = tmp ? cmatrix_multiply(O_dag, tmp) : NULL;

    cmatrix_free(tmp);

    tmp = cmatrix_multiply(raw->Sp_end, O);
    Sp_trunc = tmp ? cmatrix_multiply(O_dag, tmp) : NULL;

    cmatrix_free(tmp);
  }

  cmatrix_free(O);
  cmatrix_free(O_dag);

  if (!H_trunc || !Sz_trunc || !Sp_trunc) {
    cmatrix_free(H_trunc);
    cmatrix_free(Sz_trunc);
    cmatrix_free(Sp_trunc);

    return NULL;
  }

  dmrg_block_t *nb = malloc(sizeof(dmrg_block_t));
  if (!nb) {
    cmatrix_free(H_trunc);
    cmatrix_free(Sz_trunc);
    cmatrix_free(Sp_trunc);

    return NULL;
  }

  nb->dim = m;
  nb->H = H_trunc;
  nb->Sz_end = Sz_trunc;
  nb->Sp_end = Sp_trunc;

  return nb;
}

/*
 * One diagonalize+truncate step: enlarge cache[l-1], diagonalize
 * againstcache[r] as environment, truncate to m_max, and store into cache[l]
 * Read-only in `env`; on failure cache[l] is left as it was and 0 is returned
 */
static int sweep_step(dmrg_block_t **cache, int l, int r, double Jz, double Jxy,
                      const cmatrix_t *sz, const cmatrix_t *sp,
                      const cmatrix_t *sm, int m_max, double *E0_out,
                      double *trunc_err_out) {
  dmrg_block_t *eb = enlarge_block(cache[l - 1], Jz, Jxy, sz, sp, sm);
  if (!eb) {
    return 0;
  }

  double E0;
  cmatrix_t *Psi = lanczos_ground_state_general(eb, cache[r], Jz, Jxy, &E0);
  if (!Psi) {
    dmrg_block_free(eb);

    return 0;
  }

  double terr;
  dmrg_block_t *tb = truncate_block(eb, Psi, m_max, &terr);
  cmatrix_free(Psi);
  dmrg_block_free(eb);
  if (!tb) {
    return 0;
  }

  dmrg_block_free(cache[l]);
  cache[l] = tb;
  *E0_out = E0;
  *trunc_err_out = terr;

  return 1;
}

/*
 * Read-only version of sweep_step: same diagonalization, but truncated block is
 * discarded rather than stored, for convergence diagnostics that must not
 * perturb the cache mid-sweep
 */
static int diag_step(dmrg_block_t **cache, int l, int r, double Jz, double Jxy,
                     const cmatrix_t *sz, const cmatrix_t *sp,
                     const cmatrix_t *sm, int m_max, double *E0_out,
                     double *trunc_err_out) {
  dmrg_block_t *eb = enlarge_block(cache[l - 1], Jz, Jxy, sz, sp, sm);
  if (!eb) {
    return 0;
  }

  double E0;
  cmatrix_t *Psi = lanczos_ground_state_general(eb, cache[r], Jz, Jxy, &E0);
  if (!Psi) {
    dmrg_block_free(eb);

    return 0;
  }

  double terr;
  dmrg_block_t *tb = truncate_block(eb, Psi, m_max, &terr);
  cmatrix_free(Psi);
  dmrg_block_free(eb);
  if (!tb) {
    return 0;
  }

  dmrg_block_free(tb); // discard: diagnostic only, cache untouched
  *E0_out = E0;
  *trunc_err_out = terr;

  return 1;
}

static void free_cache(dmrg_block_t **cache, int n) {
  if (!cache) {
    return;
  }

  for (int i = 0; i < n; i++) {
    dmrg_block_free(cache[i]);
  }
  free(cache);
}

finite_dmrg_result_t *finite_dmrg_run(int N_target, double Jz, double Jxy,
                                      int m_max, int n_sweeps) {
  if (N_target < 4 || m_max < 1 || n_sweeps < 1 || !isfinite(Jz) ||
      !isfinite(Jxy)) {
    return NULL;
  }

  int N = (N_target % 2 == 0) ? N_target : N_target + 1; // round up to even
  int half = N / 2;

  cmatrix_t *sz = site_sz();
  cmatrix_t *sp = site_sp();
  cmatrix_t *sm = sp ? cmatrix_adjoint(sp) : NULL;
  dmrg_block_t **cache = calloc((size_t)N, sizeof(dmrg_block_t *));
  double *sweep_energy = malloc(sizeof(double) * (size_t)n_sweeps);
  finite_dmrg_result_t *result = malloc(sizeof(finite_dmrg_result_t));
  if (!sz || !sp || !sm || !cache || !sweep_energy || !result) {
    cmatrix_free(sz);
    cmatrix_free(sp);
    cmatrix_free(sm);
    free_cache(cache, N);
    free(sweep_energy);
    free(result);

    return NULL;
  }

  cache[1] = dmrg_block_init();
  if (!cache[1]) {
    goto fail;
  }

  // Warmup: mirrored infinite algorithm builds cache[1..half-1], each truncated
  // to m_max as usual
  int L = 1;
  while (L + 1 < half) {
    dmrg_block_t *eb = enlarge_block(cache[L], Jz, Jxy, sz, sp, sm);
    if (!eb) {
      goto fail;
    }

    double E0, terr;
    cmatrix_t *Psi = lanczos_ground_state_general(eb, eb, Jz, Jxy, &E0);
    if (!Psi) {
      dmrg_block_free(eb);

      goto fail;
    }

    dmrg_block_t *tb = truncate_block(eb, Psi, m_max, &terr);
    cmatrix_free(Psi);
    dmrg_block_free(eb);
    if (!tb) {
      goto fail;
    }

    cache[L + 1] = tb;
    L++;
  }

  // Final warmup step reaching cache[half]: also truncated to m_max
  {
    dmrg_block_t *eb = enlarge_block(cache[half - 1], Jz, Jxy, sz, sp, sm);
    if (!eb) {
      goto fail;
    }

    double E0, terr;
    cmatrix_t *Psi = lanczos_ground_state_general(eb, eb, Jz, Jxy, &E0);
    if (!Psi) {
      dmrg_block_free(eb);

      goto fail;
    }

    dmrg_block_t *tb = truncate_block(eb, Psi, m_max, &terr);
    cmatrix_free(Psi);
    dmrg_block_free(eb);
    if (!tb) {
      goto fail;
    }

    cache[half] = tb;
  }

  // NOTE: Sweep: alternately grow one side from 2 up to N-1, then other,
  // re-truncating the growing side at every step from the cache's current
  // (increasingly well-optimized) opposite-side blocks
  int direction = 1; // +1: l runs 2..N-1 upward; -1: l runs N-2..1 downward
  double last_E0 = 0.0, last_terr = 0.0;
  for (int s = 0; s < n_sweeps; s++) {
    int l_start = (direction == 1) ? 2 : N - 2;
    int l_end = (direction == 1) ? N - 1 : 1;

    for (int l = l_start; (direction == 1) ? (l <= l_end) : (l >= l_end);
         l += direction) {
      int r = N - l;
      if (r < 1 || !cache[l - 1] || !cache[r]) {
        continue;
      }

      double E0, terr;
      if (!sweep_step(cache, l, r, Jz, Jxy, sz, sp, sm, m_max, &E0, &terr)) {
        goto fail;
      }
    }

    direction = -direction;

    // Read-only diagnostic checkpoint at the most balanced bipartition
    if (!diag_step(cache, half, N - half, Jz, Jxy, sz, sp, sm, m_max, &last_E0,
                   &last_terr)) {
      goto fail;
    }
    sweep_energy[s] = last_E0;
  }

  result->energy = last_E0;
  result->energy_per_site = last_E0 / (double)N;
  result->N_reached = N;
  result->m_max = m_max;
  result->n_sweeps = n_sweeps;
  result->truncation_error = last_terr;
  result->sweep_energy = sweep_energy;

  cmatrix_free(sz);
  cmatrix_free(sp);
  cmatrix_free(sm);
  free_cache(cache, N);

  return result;

fail:
  cmatrix_free(sz);
  cmatrix_free(sp);
  cmatrix_free(sm);
  free_cache(cache, N);
  free(sweep_energy);
  free(result);

  return NULL;
}

void finite_dmrg_result_free(finite_dmrg_result_t *r) {
  if (!r) {
    return;
  }

  free(r->sweep_energy);
  free(r);
}
