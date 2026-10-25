/*
 * Coupled Cluster Singles and Doubles (CCSD)
 *
 * Equation-of-motion ansatz (Reference: Purvis & Bartlett 1982; Gauss & Stanton
 * in "Advances in Chemical Physics" 123 for diagrammatic formulation used
 * below):
 *
 *   |\Psi_CCSD> = \exp(T1 + T2) |\Phi_0>
 *
 * Where
 *  |\Phi_0> is closed-shell RHF reference and
 *
 *   T1 = \sum_{i,a} t_i^a c^\dagger_a c_i
 *   T2 = 1/4 \sum_{i,j,a,b} t_{ij}^{ab} c^\dagger_a c^\dagger_b c_j c_i
 *
 * are cluster operators in second quantization
 * Correlation energy is
 *
 *   E_corr = <\Phi_0| (H_N \exp(T1 + T2))_c |\Phi_0>
 *
 * and amplitudes solve projected equations
 *
 *   <\Phi_i^a| H_N \exp(T1 + T2) |\Phi_0>_c = 0
 *   <\Phi_{ij}^{ab}| H_N \exp(T1 + T2) |\Phi_0>_c = 0
 *
 * iterated to convergence with a fixed-point (DIIS-free) update. Since
 * CCSD is exact through fourth order in MP perturbation theory and its
 * energy is a sum of terms up to (T2)^2, amplitude equations terminate
 * at finite diagrams
 *
 * Everything is in spin-orbital basis (nso = 2 * n_spatial):
 *   - MO integrals come in as spatial chemist-notation (pq|rs) on n_spatial
 *     orbitals; spin-orbital physicist-notation antisymmetrized tensor
 *     <pq||rs> = (pr|qs) - (ps|qr) is evaluated on demand via v_elem, so
 *     nothing nso^4 = 16 n_spatial^4 ever has to be resident during iterative
 *     solve
 *   - occ[] and virt[] are explicit spin-orbital index lists over active window
 *     (frozen-core orbitals excluded from both)
 *   - Fso is diagonal (canonical RHF), so amplitude energy denominator
 *     d_ij^ab = F_ii + F_jj - F_aa - F_bb is built from 4 scalars
 *   - 2 T1-dependent τ intermediates are antisymmetrized products
 *     \tau_{ij}^{ab}         = t_{ij}^{ab} + t_i^a t_j^b - t_i^b t_j^a
 *     \tilde{\tau}_{ij}^{ab} = t_{ij}^{ab} + 1/2 (t_i^a t_j^b - t_i^b t_j^a)
 *     used in standard CCSD diagrams
 *
 * Intermediates: Each iteration rebuilds Fae, Fmi, Fme, Wmnij and Wmbej from
 * current (T1, T2); Wabef is never stored, only computed on demand for (ab|ef)
 * contractions of T2 equation. This mirrors standard "5-intermediate" CCSD
 * factorization and keeps persistent memory budget at O(nso^2) plus T2 tensor
 * itself, at cost of recomputing a few of smaller intermediates once per outer
 * RHS evaluation
 *
 * Convergence: MP2 amplitude is used as initial guess. Iteration is declared
 * converged when correlation energy changes by less than conv_tol between two
 * consecutive sweeps (no amplitude-norm criterion). At second order MP2 alone
 * would already be an upper bound of CCSD correlation energy; outer loop
 * reaches CCSD fixed point in a handful of iterations for small systems this
 * module targets (typically H2/STO-3G through LiH/STO-3G, n_spatial <= 8)
 *
 * Amplitudes exit: ccsd_run_ex() can optionally hand back converged amplitudes
 * for post-CCSD methods like perturbative triples correction in ccsd_t.c. Exit
 * path materializes full antisymmetrized V tensor here (it is only nso^4 for
 * small active spaces this is used on) so that perturbative_triples() can index
 * it directly without paying v_elem's cost inside innermost triple loop
 */

#include "ccsd.h"
#include "molecular_integrals.h"
#include <math.h>
#include <stdlib.h>

// Flat 2-index / 4-index accessors over spin-orbital dimension nso
#define IDX2(arr, nso, p, q) (arr)[(size_t)(p) * (nso) + (q)]
#define IDX4(arr, nso, p, q, r, s)                                             \
  (arr)[(((size_t)(p) * (nso) + (q)) * (nso) + (r)) * (nso) + (s)]

static int spin_of(int orb_idx) { return orb_idx % 2; }

static int spatial_of(int orb_idx) { return orb_idx / 2; }

// Single-element on-the-fly version of build_antisymmetrized (below): same
// <pq||rs> formula, computed directly from n_spatial-sized chemist eri_mo
// tensor for one (p,q,r,s)
//
// NOTE: Antisymmetrized integrals never need nso^4 = 16 * n_spatial^4 of
// persistent storage during iterative solve - only n_spatial^4 is ever
// resident
static inline double v_elem(int n_spatial, const double *eri_mo, int p_idx,
                            int q_idx, int r_idx, int s_idx) {
  int spin_p = spin_of(p_idx), spatial_p = spatial_of(p_idx);
  int spin_q = spin_of(q_idx), spatial_q = spatial_of(q_idx);
  int spin_r = spin_of(r_idx), spatial_r = spatial_of(r_idx);
  int spin_s = spin_of(s_idx), spatial_s = spatial_of(s_idx);

  double v1_val = (spin_p == spin_r && spin_q == spin_s)
                      ? MOLINT_ERI(eri_mo, n_spatial, spatial_p, spatial_r,
                                   spatial_q, spatial_s)
                      : 0.0;
  double v2_val = (spin_p == spin_s && spin_q == spin_r)
                      ? MOLINT_ERI(eri_mo, n_spatial, spatial_p, spatial_s,
                                   spatial_q, spatial_r)
                      : 0.0;

  return v1_val - v2_val;
}

// Antisymmetrized physicist-notation spin-orbital integrals from spatial
// chemist-notation MO integrals:
//   <pq|rs> = (pr|qs)                    [chemist -> physicist]
//   <pq||rs> = <pq|rs> - <pq|sr>         [antisymmetrize]
// with spin conservation (pq|rs)'s spin-orbital version is nonzero only when
// spin_p = spin_r and spin_q = spin_s (for <pq|rs>), or spin_p = spin_s and
// spin_q = spin_r (for <pq|sr>)
static double *build_antisymmetrized(int n_spatial, const double *eri_mo,
                                     int nso) {
  double *v_tensor =
      (double *)calloc((size_t)nso * nso * nso * nso, sizeof(double));
  if (!v_tensor) {
    return NULL;
  }

  for (int p_idx = 0; p_idx < nso; p_idx++) {
    int spin_p = spin_of(p_idx);
    int spatial_p = spatial_of(p_idx);

    for (int q_idx = 0; q_idx < nso; q_idx++) {
      int spin_q = spin_of(q_idx);
      int spatial_q = spatial_of(q_idx);

      for (int r_idx = 0; r_idx < nso; r_idx++) {
        int spin_r = spin_of(r_idx);
        int spatial_r = spatial_of(r_idx);

        for (int s_idx = 0; s_idx < nso; s_idx++) {
          int spin_s = spin_of(s_idx);
          int spatial_s = spatial_of(s_idx);

          double v1_val = (spin_p == spin_r && spin_q == spin_s)
                              ? MOLINT_ERI(eri_mo, n_spatial, spatial_p,
                                           spatial_r, spatial_q, spatial_s)
                              : 0.0;
          double v2_val = (spin_p == spin_s && spin_q == spin_r)
                              ? MOLINT_ERI(eri_mo, n_spatial, spatial_p,
                                           spatial_s, spatial_q, spatial_r)
                              : 0.0;

          IDX4(v_tensor, nso, p_idx, q_idx, r_idx, s_idx) = v1_val - v2_val;
        }
      }
    }
  }

  return v_tensor;
}

// \tau_ij^ab = t2_ijab + t1_ia * t1_jb - t1_ib * t1_ja
// (full antisymmetrized T1xT1 product, used in Wmnij/Wabef/T2-update)
static double tau_full(const double *t1_amp, const double *t2_amp, int nso,
                       int occ_i, int occ_j, int virt_a, int virt_b) {
  return IDX4(t2_amp, nso, occ_i, occ_j, virt_a, virt_b) +
         IDX2(t1_amp, nso, occ_i, virt_a) * IDX2(t1_amp, nso, occ_j, virt_b) -
         IDX2(t1_amp, nso, occ_i, virt_b) * IDX2(t1_amp, nso, occ_j, virt_a);
}

// \tau-\tilde_ij^ab = t2_ijab + 0.5 * (t1_ia * t1_jb - t1_ib * t1_ja)
static double tau_tilde(const double *t1_amp, const double *t2_amp, int nso,
                        int occ_i, int occ_j, int virt_a, int virt_b) {
  return IDX4(t2_amp, nso, occ_i, occ_j, virt_a, virt_b) +
         0.5 * (IDX2(t1_amp, nso, occ_i, virt_a) *
                    IDX2(t1_amp, nso, occ_j, virt_b) -
                IDX2(t1_amp, nso, occ_i, virt_b) *
                    IDX2(t1_amp, nso, occ_j, virt_a));
}

typedef struct {
  int nso;
  int nocc, nvirt;
  int *occ, *virt;      /* nocc / nvirt long, spin-orbital indices */
  double *Fso;          /* nso x nso, diagonal (canonical orbitals) */
  double *Dia;          /* nso x nso, only occ x virt entries used */
  int n_spatial;        /* for v_elem(): number of spatial MOs */
  const double *eri_mo; /* for v_elem(): caller-owned, not freed here */
} ccsd_ctx_t;

// d_ijab(i,j,a,b) = Fso[i] + Fso[j] - Fso[a] - Fso[b]
// all diagonal, canonical orbitals - computed on demand rather than cached in
// an nso^4 array
static inline double d_ijab(const ccsd_ctx_t *ctx, int i, int j, int a, int b) {
  int nso = ctx->nso;
  const double *Fso = ctx->Fso;

  return IDX2(Fso, nso, i, i) + IDX2(Fso, nso, j, j) - IDX2(Fso, nso, a, a) -
         IDX2(Fso, nso, b, b);
}

// Build standard CCSD one- and two-body intermediates from current amplitudes
// (t1, t2) in spin-orbital basis
//
// This routine assembles intermediates used by CCSD amplitude updates:
// - Fae(a,e)      : dressed virtual-virtual Fock block
// - Fmi(m,i)      : dressed occupied-occupied Fock block
// - Fme(m,e)      : occupied-virtual coupling block
// - Wmnij(m,n,i,j): occupied-occupied two-body intermediate
// - Wmbej(m,b,e,j): mixed occupied/virtual two-body intermediate
//
// Inputs:
// - ctx: orbital partitions (occ/virt), spin-orbital Fock matrix, and MO ERIs
// - t1, t2: current singles/doubles amplitudes
//
// Outputs (caller-allocated):
// - Fae, Fmi, Fme, Wmnij, Wmbej (all overwritten by this routine)
//
// Indexing / conventions:
// - occ[] and virt[] store explicit spin-orbital indices
// - i,j,m,n refer to occupied spin orbitals; a,b,e,f refer to virtual ones
// - V(p,q,r,s) expands to <pq||rs> evaluated on demand from spatial ERIs with
//   spin selection rules (via v_elem), avoiding persistent nso^4 storage
//
// NOTE:: Routine is intentionally monolithic for performance (tight nested
// loops and minimal temporary allocations)
static void build_intermediates(const ccsd_ctx_t *ctx, const double *t1,
                                const double *t2, double *Fae, double *Fmi,
                                double *Fme, double *Wmnij, double *Wmbej) {
  int nso = ctx->nso;
  const int *occ = ctx->occ, *virt = ctx->virt;
  int no = ctx->nocc, nv = ctx->nvirt;
  const double *Fso = ctx->Fso;
  int n_spatial = ctx->n_spatial;
  const double *eri_mo = ctx->eri_mo;

#define V(p, q, r, s) v_elem(n_spatial, eri_mo, (p), (q), (r), (s))

// Build Fae(a,e): virtual-virtual dressed Fock-like intermediate
// Structure mirrors standard CCSD equations:
//   - linear singles/Fock coupling
//   - singles-driven two-electron dressing
//   - \tau_\tilde (T2 + disconnected T1T1) correction
#pragma omp parallel for schedule(dynamic)
  for (int ai = 0; ai < nv; ai++) {
    int a = virt[ai];

    for (int ei = 0; ei < nv; ei++) {
      int e = virt[ei];
      double s = 0.0;

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];
        // -1/2 * f_me * t_m^a contribution
        s -= 0.5 * IDX2(Fso, nso, m, e) * IDX2(t1, nso, m, a);

        for (int fi = 0; fi < nv; fi++) {
          int f = virt[fi];
          // + t_m^f * <ma||fe>
          s += IDX2(t1, nso, m, f) * V(m, a, f, e);
        }

        for (int ni = 0; ni < no; ni++) {
          int n = occ[ni];

          for (int fi = 0; fi < nv; fi++) {
            int f = virt[fi];
            // -1/2 * tau_tilde_mn^{af} * <mn||ef>
            s -= 0.5 * tau_tilde(t1, t2, nso, m, n, a, f) * V(m, n, e, f);
          }
        }
      }

      // Diagonal uses 0.0 here (canonical convention in this implementation),
      // then add accumulated dressing terms.
      IDX2(Fae, nso, a, e) = (a != e ? IDX2(Fso, nso, a, e) : 0.0) + s;
    }
  }

// Build Fmi(m,i): occupied-occupied dressed Fock-like intermediate
// Structure mirrors standard CCSD equations:
//   - linear singles/Fock coupling
//   - singles-driven two-electron dressing
//   - \tau_\tilde (T2 + disconnected T1T1) correction
#pragma omp parallel for schedule(dynamic)
  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    for (int ii = 0; ii < no; ii++) {
      int i = occ[ii];
      double s = 0.0;

      for (int ei = 0; ei < nv; ei++) {
        int e = virt[ei];
        // +1/2 * f_me * t_i^e contribution
        s += 0.5 * IDX2(t1, nso, i, e) * IDX2(Fso, nso, m, e);

        for (int ni = 0; ni < no; ni++) {
          int n = occ[ni];
          // + t_n^e * <mn||ie>
          s += IDX2(t1, nso, n, e) * V(m, n, i, e);
        }

        for (int ni = 0; ni < no; ni++) {
          int n = occ[ni];

          for (int fi = 0; fi < nv; fi++) {
            int f = virt[fi];
            // +1/2 * tau_tilde_in^{ef} * <mn||ef>
            s += 0.5 * tau_tilde(t1, t2, nso, i, n, e, f) * V(m, n, e, f);
          }
        }
      }

      IDX2(Fmi, nso, m, i) = (m != i ? IDX2(Fso, nso, m, i) : 0.0) + s;
    }
  }

// Build Fme(m,e): occupied-virtual coupling intermediate
// Structure:
//   - bare Fock matrix element f_me
//   - T1-driven two-electron dressing <mn||ef>
#pragma omp parallel for schedule(dynamic)
  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    for (int ei = 0; ei < nv; ei++) {
      int e = virt[ei];
      // Bare Fock element f_me
      double s = IDX2(Fso, nso, m, e);

      for (int ni = 0; ni < no; ni++) {
        int n = occ[ni];

        for (int fi = 0; fi < nv; fi++) {
          int f = virt[fi];
          // + t_n^f * <mn||ef>
          s += IDX2(t1, nso, n, f) * V(m, n, e, f);
        }
      }

      IDX2(Fme, nso, m, e) = s;
    }
  }

// Build Wmnij(m,n,i,j): occupied-occupied two-body intermediate
// Structure:
//   - bare antisymmetrized integral <mn||ij>
//   - linear singles dressing
//   - \tau (T2 + full T1T1) dressing
#pragma omp parallel for schedule(dynamic)
  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    for (int ni = 0; ni < no; ni++) {
      int n = occ[ni];

      for (int ii = 0; ii < no; ii++) {
        int i = occ[ii];

        for (int ji = 0; ji < no; ji++) {
          int j = occ[ji];
          // Bare integral <mn||ij>
          double s = V(m, n, i, j);

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];
            // + t_j^e <mn||ie> - t_i^e <mn||je>
            s += IDX2(t1, nso, j, e) * V(m, n, i, e) -
                 IDX2(t1, nso, i, e) * V(m, n, j, e);
          }

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];

            for (int fi = 0; fi < nv; fi++) {
              int f = virt[fi];
              // + 1/4 * tau_full_ij^{ef} * <mn||ef>
              s += 0.25 * tau_full(t1, t2, nso, i, j, e, f) * V(m, n, e, f);
            }
          }

          IDX4(Wmnij, nso, m, n, i, j) = s;
        }
      }
    }
  }

// Build Wmbej(m,b,e,j): mixed occupied/virtual two-body intermediate
// Structure:
//   - bare antisymmetrized integral <mb||ej>
//   - linear singles dressing (both virtual and occupied)
//   - quadratic T2 + T1T1 dressing
#pragma omp parallel for schedule(dynamic)
  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    for (int bi = 0; bi < nv; bi++) {
      int b = virt[bi];

      for (int ei = 0; ei < nv; ei++) {
        int e = virt[ei];

        for (int ji = 0; ji < no; ji++) {
          int j = occ[ji];
          // Bare integral <mb||ej>
          double s = V(m, b, e, j);

          for (int fi = 0; fi < nv; fi++) {
            int f = virt[fi];
            // + t_j^f * <mb||ef>
            s += IDX2(t1, nso, j, f) * V(m, b, e, f);
          }

          for (int ni = 0; ni < no; ni++) {
            int n = occ[ni];
            // - t_n^b * <mn||ej>
            s -= IDX2(t1, nso, n, b) * V(m, n, e, j);
          }

          for (int ni = 0; ni < no; ni++) {
            int n = occ[ni];

            for (int fi = 0; fi < nv; fi++) {
              int f = virt[fi];
              // - (1/2 t_jn^{fb} + t_j^f t_n^b) * <mn||ef>
              s -= (0.5 * IDX4(t2, nso, j, n, f, b) +
                    IDX2(t1, nso, j, f) * IDX2(t1, nso, n, b)) *
                   V(m, n, e, f);
            }
          }

          IDX4(Wmbej, nso, m, b, e, j) = s;
        }
      }
    }
  }

#undef V
}

// Single element of Wabef, computed on demand rather than cached in a
// persistent nv^4 array. Wabef-building loop body, evaluated for one (a,b,e,f)
// at point of use
static double wabef_elem(const ccsd_ctx_t *ctx, const double *t1,
                         const double *t2, int a, int b, int e, int f) {
  int nso = ctx->nso;
  const int *occ = ctx->occ;
  int no = ctx->nocc;
  int n_spatial = ctx->n_spatial;
  const double *eri_mo = ctx->eri_mo;

#define V(p, q, r, s) v_elem(n_spatial, eri_mo, (p), (q), (r), (s))

  // Bare integral <ab||ef>
  double s = V(a, b, e, f);

  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    // - t_m^b <am||ef> + t_m^a <bm||ef>
    s -= IDX2(t1, nso, m, b) * V(a, m, e, f) -
         IDX2(t1, nso, m, a) * V(b, m, e, f);
  }

  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    for (int ni = 0; ni < no; ni++) {
      int n = occ[ni];
      // + 1/4 * tau_full_ab^{mn} * <mn||ef>
      s += 0.25 * tau_full(t1, t2, nso, m, n, a, b) * V(m, n, e, f);
    }
  }

#undef V

  return s;
}

// Perform one CCSD amplitude update step
//
// Inputs:
// - ctx: system dimensions, orbital partitions, Fock matrix, and ERIs
// - t1, t2: current singles/doubles amplitudes
//
// Outputs:
// - new_t1, new_t2: next-iteration amplitudes (same indexing as t1/t2)
// - Fae, Fmi, Fme, Wmnij, Wmbej: reusable intermediates populated for this step
//
// High-level flow:
// 1. Build CCSD intermediates from current amplitudes
// 2. Update T1 amplitudes over occupied/virtual pairs
// 3. Update T2 amplitudes over occupied/occupied/virtual/virtual tuples
// 4. Apply orbital-energy denominators to form final residual updates
//
// NOTE: Heavy tensor contractions are parallelized with OpenMP
// Kept antisymmetrized two-electron integrals on-demand through V(...) to avoid
// persistent nso^4 storage in iterative path
static void update_amplitudes(const ccsd_ctx_t *ctx, const double *t1,
                              const double *t2, double *new_t1, double *new_t2,
                              double *Fae, double *Fmi, double *Fme,
                              double *Wmnij, double *Wmbej) {
  int nso = ctx->nso;
  const int *occ = ctx->occ, *virt = ctx->virt;
  int no = ctx->nocc, nv = ctx->nvirt;
  const double *Fso = ctx->Fso;
  int n_spatial = ctx->n_spatial;
  const double *eri_mo = ctx->eri_mo;

#define V(p, q, r, s) v_elem(n_spatial, eri_mo, (p), (q), (r), (s))

  // Build reusable intermediates once per iteration
  build_intermediates(ctx, t1, t2, Fae, Fmi, Fme, Wmnij, Wmbej);

// T1 update:
// - assemble T1 residual for each (i,a)
// - divide by energy denominator at end
// Structure mirrors standard CCSD T1 diagrams:
//   - bare f_ia
//   - Fock dressings with singles (Fae, Fmi)
//   - T2 x Fme coupling
//   - bare <na||if> term
//   - T2 x bare integral ring terms
#pragma omp parallel for schedule(dynamic)
  for (int ii = 0; ii < no; ii++) {
    int i = occ[ii];

    for (int ai = 0; ai < nv; ai++) {
      int a = virt[ai];
      // Bare f_ia
      double s = IDX2(Fso, nso, i, a);

      for (int ei = 0; ei < nv; ei++) {
        int e = virt[ei];
        // + t_i^e * Fae_ae
        s += IDX2(t1, nso, i, e) * IDX2(Fae, nso, a, e);
      }

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];
        // - t_m^a * Fmi_mi
        s -= IDX2(t1, nso, m, a) * IDX2(Fmi, nso, m, i);
      }

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];

        for (int ei = 0; ei < nv; ei++) {
          int e = virt[ei];
          // + t_im^{ae} * Fme_me
          s += IDX4(t2, nso, i, m, a, e) * IDX2(Fme, nso, m, e);
        }
      }

      for (int ni = 0; ni < no; ni++) {
        int n = occ[ni];

        for (int fi = 0; fi < nv; fi++) {
          int f = virt[fi];
          // - t_n^f * <na||if>
          s -= IDX2(t1, nso, n, f) * V(n, a, i, f);
        }
      }

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];

        for (int ei = 0; ei < nv; ei++) {
          int e = virt[ei];

          for (int fi = 0; fi < nv; fi++) {
            int f = virt[fi];
            // - 1/2 * t_im^{ef} * <ma||ef>
            s -= 0.5 * IDX4(t2, nso, i, m, e, f) * V(m, a, e, f);
          }
        }
      }

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];

        for (int ei = 0; ei < nv; ei++) {
          int e = virt[ei];

          for (int ni = 0; ni < no; ni++) {
            int n = occ[ni];
            // - 1/2 * t_mn^{ae} * <nm||ei>
            s -= 0.5 * IDX4(t2, nso, m, n, a, e) * V(n, m, e, i);
          }
        }
      }

      IDX2(new_t1, nso, i, a) = s / IDX2(ctx->Dia, nso, i, a);
    }
  }

// T2 update:
// - assemble T2 residual for each (i,j,a,b)
// - divide by energy denominator at end
// Structure mirrors standard CCSD T2 diagrams:
//   - bare <ij||ab>
//   - Fock dressings with T2 (Fae, Fmi)
//   - \tau x Wmnij / Wabef ladder terms
//   - Wmbej crossed-ring term
//   - T1 x bare integral linear terms
#pragma omp parallel for schedule(dynamic)
  for (int ii = 0; ii < no; ii++) {
    int i = occ[ii];

    for (int ji = 0; ji < no; ji++) {
      int j = occ[ji];

      for (int ai = 0; ai < nv; ai++) {
        int a = virt[ai];

        for (int bi = 0; bi < nv; bi++) {
          int b = virt[bi];
          // Bare integral <ij||ab>
          double s = V(i, j, a, b);

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];
            double corr_b = 0.0;
            double corr_a = 0.0;

            for (int mi = 0; mi < no; mi++) {
              int m = occ[mi];
              corr_b += IDX2(t1, nso, m, b) * IDX2(Fme, nso, m, e);
              corr_a += IDX2(t1, nso, m, a) * IDX2(Fme, nso, m, e);
            }

            // + t_ij^{ae} * (Fae_be - 1/2 t_m^b Fme_me)
            // - t_ij^{be} * (Fae_ae - 1/2 t_m^a Fme_me)
            s += IDX4(t2, nso, i, j, a, e) *
                     (IDX2(Fae, nso, b, e) - 0.5 * corr_b) -
                 IDX4(t2, nso, i, j, b, e) *
                     (IDX2(Fae, nso, a, e) - 0.5 * corr_a);
          }

          for (int mi = 0; mi < no; mi++) {
            int m = occ[mi];
            double corr_j = 0.0, corr_i = 0.0;

            for (int ei = 0; ei < nv; ei++) {
              int e = virt[ei];
              corr_j += IDX2(t1, nso, j, e) * IDX2(Fme, nso, m, e);
              corr_i += IDX2(t1, nso, i, e) * IDX2(Fme, nso, m, e);
            }

            // - t_im^{ab} * (Fmi_mj + 1/2 t_j^e Fme_me)
            // + t_jm^{ab} * (Fmi_mi + 1/2 t_i^e Fme_me)
            s -= IDX4(t2, nso, i, m, a, b) *
                     (IDX2(Fmi, nso, m, j) + 0.5 * corr_j) -
                 IDX4(t2, nso, j, m, a, b) *
                     (IDX2(Fmi, nso, m, i) + 0.5 * corr_i);
          }

          for (int mi = 0; mi < no; mi++) {
            int m = occ[mi];

            for (int ni = 0; ni < no; ni++) {
              int n = occ[ni];
              // + 1/2 * tau_mn^{ab} * Wmnij_mnij
              s += 0.5 * tau_full(t1, t2, nso, m, n, a, b) *
                   IDX4(Wmnij, nso, m, n, i, j);
            }
          }

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];

            for (int fi = 0; fi < nv; fi++) {
              int f = virt[fi];
              // + 1/2 * tau_ij^{ef} * Wabef_abef
              s += 0.5 * tau_full(t1, t2, nso, i, j, e, f) *
                   wabef_elem(ctx, t1, t2, a, b, e, f);
            }
          }

          for (int mi = 0; mi < no; mi++) {
            int m = occ[mi];

            for (int ei = 0; ei < nv; ei++) {
              int e = virt[ei];
              // Crossed-ring Wmbej terms with T2 and T1T1:
              // + t_im^{ae} Wmbej_mbej - t_i^e t_m^a <mb||ej>
              // - t_im^{be} Wmbej_maej + t_i^e t_m^b <ma||ej>
              // - t_jm^{ae} Wmbej_mbei + t_j^e t_m^a <mb||ei>
              // + t_jm^{be} Wmbej_maei - t_j^e t_m^b <ma||ei>
              s += IDX4(t2, nso, i, m, a, e) * IDX4(Wmbej, nso, m, b, e, j) -
                   IDX2(t1, nso, i, e) * IDX2(t1, nso, m, a) * V(m, b, e, j);
              s -= IDX4(t2, nso, i, m, b, e) * IDX4(Wmbej, nso, m, a, e, j) -
                   IDX2(t1, nso, i, e) * IDX2(t1, nso, m, b) * V(m, a, e, j);
              s -= IDX4(t2, nso, j, m, a, e) * IDX4(Wmbej, nso, m, b, e, i) -
                   IDX2(t1, nso, j, e) * IDX2(t1, nso, m, a) * V(m, b, e, i);
              s += IDX4(t2, nso, j, m, b, e) * IDX4(Wmbej, nso, m, a, e, i) -
                   IDX2(t1, nso, j, e) * IDX2(t1, nso, m, b) * V(m, a, e, i);
            }
          }

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];
            // + t_i^e <ab||ej> - t_j^e <ab||ei>
            s += IDX2(t1, nso, i, e) * V(a, b, e, j) -
                 IDX2(t1, nso, j, e) * V(a, b, e, i);
          }

          for (int mi = 0; mi < no; mi++) {
            int m = occ[mi];
            // - t_m^a <mb||ij> + t_m^b <ma||ij>
            s -= IDX2(t1, nso, m, a) * V(m, b, i, j) -
                 IDX2(t1, nso, m, b) * V(m, a, i, j);
          }

          IDX4(new_t2, nso, i, j, a, b) = s / d_ijab(ctx, i, j, a, b);
        }
      }
    }
  }

#undef V
}

// Compute CCSD correlation energy:
//   E_corr = \sum_{ia} f_ia t_i^a
//          + 1/4 \sum_{ijab} <ij||ab> t_ij^{ab}
//          + 1/2 \sum_{ijab} <ij||ab> t_i^a t_j^b
// This is standard CCSD energy expression (H_N * T1 * T2 diagram contributes
// last term)
static double ccsd_energy(const ccsd_ctx_t *ctx, const double *t1,
                          const double *t2) {
  int nso = ctx->nso;
  const int *occ = ctx->occ, *virt = ctx->virt;
  int no = ctx->nocc, nv = ctx->nvirt;
  const double *Fso = ctx->Fso;
  int n_spatial = ctx->n_spatial;
  const double *eri_mo = ctx->eri_mo;

#define V(p, q, r, s) v_elem(n_spatial, eri_mo, (p), (q), (r), (s))

  double E = 0.0;
  for (int ii = 0; ii < no; ii++) {
    int i = occ[ii];

    for (int ai = 0; ai < nv; ai++) {
      int a = virt[ai];
      // + f_ia t_i^a (Brillouin's theorem makes this zero for canonical RHF,
      // but kept for generality)
      E += IDX2(Fso, nso, i, a) * IDX2(t1, nso, i, a);
    }
  }

  for (int ii = 0; ii < no; ii++) {
    int i = occ[ii];

    for (int ji = 0; ji < no; ji++) {
      int j = occ[ji];

      for (int ai = 0; ai < nv; ai++) {
        int a = virt[ai];

        for (int bi = 0; bi < nv; bi++) {
          int b = virt[bi];
          double v = V(i, j, a, b);

          // + 1/4 <ij||ab> t_ij^{ab}
          E += 0.25 * v * IDX4(t2, nso, i, j, a, b);
          // + 1/2 <ij||ab> t_i^a t_j^b
          E += 0.5 * v * IDX2(t1, nso, i, a) * IDX2(t1, nso, j, b);
        }
      }
    }
  }

#undef V

  return E;
}

void ccsd_amplitudes_free(ccsd_amplitudes_t *amp) {
  if (!amp) {
    return;
  }

  free(amp->occ);
  free(amp->virt);
  free(amp->Fso);
  free(amp->V);
  free(amp->t1);
  free(amp->t2);
  free(amp);
}

ccsd_result_t *ccsd_run(int n_spatial, const double *h_mo, const double *eri_mo,
                        const double *mo_energy, int n_electrons,
                        int n_frozen_spatial, double e_rhf, double conv_tol,
                        int max_iter) {
  return ccsd_run_ex(n_spatial, h_mo, eri_mo, mo_energy, n_electrons,
                     n_frozen_spatial, e_rhf, conv_tol, max_iter, NULL);
}

ccsd_result_t *ccsd_run_ex(int n_spatial, const double *h_mo,
                           const double *eri_mo, const double *mo_energy,
                           int n_electrons, int n_frozen_spatial, double e_rhf,
                           double conv_tol, int max_iter,
                           ccsd_amplitudes_t **amplitudes_out) {
  if (amplitudes_out) {
    *amplitudes_out = NULL;
  }

  if (n_spatial <= 0 || !h_mo || !eri_mo || !mo_energy ||
      n_electrons % 2 != 0 || n_frozen_spatial < 0 ||
      2 * n_frozen_spatial >= n_electrons || n_electrons > 2 * n_spatial) {
    return NULL;
  }

  // canonical RHF: off-diagonal core Hamiltonian in MO basis is zero by
  // Brillouin's theorem, only orbital energies (diagonal) are needed
  (void)h_mo;

  int nso = 2 * n_spatial;
  ccsd_ctx_t ctx;

  ctx.nso = nso;
  ctx.nocc = n_electrons - 2 * n_frozen_spatial;
  ctx.nvirt = nso - n_electrons;
  ctx.occ = calloc((size_t)ctx.nocc, sizeof(int));
  ctx.virt = calloc((size_t)ctx.nvirt, sizeof(int));
  ctx.Fso = calloc((size_t)nso * nso, sizeof(double));
  ctx.Dia = calloc((size_t)nso * nso, sizeof(double));
  ctx.n_spatial = n_spatial;
  ctx.eri_mo = eri_mo; // caller-owned; ctx never frees this

  if (!ctx.occ || !ctx.virt || !ctx.Fso || !ctx.Dia) {
    free(ctx.occ);
    free(ctx.virt);
    free(ctx.Fso);
    free(ctx.Dia);

    return NULL;
  }

  for (int p = 0; p < nso; p++) {
    IDX2(ctx.Fso, nso, p, p) = mo_energy[spatial_of(p)];
  }

  for (int i = 0, k = 0; i < nso; i++) {
    if (i >= 2 * n_frozen_spatial && i < n_electrons) {
      ctx.occ[k++] = i;
    }
  }

  for (int a = 0, k = 0; a < nso; a++) {
    if (a >= n_electrons) {
      ctx.virt[k++] = a;
    }
  }

  int no = ctx.nocc, nv = ctx.nvirt;
  for (int ii = 0; ii < no; ii++) {
    int i = ctx.occ[ii];

    for (int ai = 0; ai < nv; ai++) {
      int a = ctx.virt[ai];

      IDX2(ctx.Dia, nso, i, a) =
          IDX2(ctx.Fso, nso, i, i) - IDX2(ctx.Fso, nso, a, a);
    }
  }

  double *t1 = calloc((size_t)nso * nso, sizeof(double));
  double *t2 = calloc((size_t)nso * nso * nso * nso, sizeof(double));
  double *new_t1 = calloc((size_t)nso * nso, sizeof(double));
  double *new_t2 = calloc((size_t)nso * nso * nso * nso, sizeof(double));
  double *Fae = calloc((size_t)nso * nso, sizeof(double));
  double *Fmi = calloc((size_t)nso * nso, sizeof(double));
  double *Fme = calloc((size_t)nso * nso, sizeof(double));
  double *Wmnij = calloc((size_t)nso * nso * nso * nso, sizeof(double));
  double *Wmbej = calloc((size_t)nso * nso * nso * nso, sizeof(double));

  ccsd_result_t *result = malloc(sizeof(ccsd_result_t));

  if (!t1 || !t2 || !new_t1 || !new_t2 || !Fae || !Fmi || !Fme || !Wmnij ||
      !Wmbej || !result) {
    free(t1);
    free(t2);
    free(new_t1);
    free(new_t2);
    free(Fae);
    free(Fmi);
    free(Fme);
    free(Wmnij);
    free(Wmbej);
    free(result);
    free(ctx.occ);
    free(ctx.virt);
    free(ctx.Fso);
    free(ctx.Dia);

    return NULL;
  }

  // MP2 initial guess
  for (int ii = 0; ii < no; ii++) {
    int i = ctx.occ[ii];

    for (int ji = 0; ji < no; ji++) {
      int j = ctx.occ[ji];

      for (int ai = 0; ai < nv; ai++) {
        int a = ctx.virt[ai];

        for (int bi = 0; bi < nv; bi++) {
          int b = ctx.virt[bi];

          IDX4(t2, nso, i, j, a, b) =
              v_elem(n_spatial, eri_mo, i, j, a, b) / d_ijab(&ctx, i, j, a, b);
        }
      }
    }
  }

  double e_prev = 0.0;
  int converged = 0, it;

  for (it = 0; it < max_iter; it++) {
    update_amplitudes(&ctx, t1, t2, new_t1, new_t2, Fae, Fmi, Fme, Wmnij,
                      Wmbej);

    double *tmp;

    tmp = t1;
    t1 = new_t1;
    new_t1 = tmp;
    tmp = t2;
    t2 = new_t2;
    new_t2 = tmp;

    double e_corr = ccsd_energy(&ctx, t1, t2);
    if (fabs(e_corr - e_prev) < conv_tol) {
      converged = 1;
      e_prev = e_corr;
      it++;

      break;
    }

    e_prev = e_corr;
  }

  result->correlation_energy = e_prev;
  result->total_energy = e_rhf + e_prev;
  result->converged = converged;
  result->iterations = it;

  if (amplitudes_out) {
    ccsd_amplitudes_t *amp = malloc(sizeof(ccsd_amplitudes_t));
    double *v_materialized =
        amp ? build_antisymmetrized(n_spatial, eri_mo, nso) : NULL;
    if (amp && v_materialized) {
      amp->nso = nso;
      amp->nocc = no;
      amp->nvirt = nv;
      amp->occ = ctx.occ;
      amp->virt = ctx.virt;
      amp->Fso = ctx.Fso;
      amp->V = v_materialized;
      amp->t1 = t1;
      amp->t2 = t2;
      *amplitudes_out = amp;
    } else {
      // allocation failure: fall back to freeing everything below
      free(amp);
      free(v_materialized);
      free(ctx.occ);
      free(ctx.virt);
      free(ctx.Fso);
      free(t1);
      free(t2);
    }
  } else {
    free(ctx.occ);
    free(ctx.virt);
    free(ctx.Fso);
    free(t1);
    free(t2);
  }

  free(new_t1);
  free(new_t2);
  free(Fae);
  free(Fmi);
  free(Fme);
  free(Wmnij);
  free(Wmbej);
  free(ctx.Dia);

  return result;
}
