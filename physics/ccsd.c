#include "ccsd.h"
#include "molecular_integrals.h"
#include <math.h>
#include <stdlib.h>

// Flat 2-index / 4-index accessors over the spin-orbital dimension nso
#define IDX2(arr, nso, p, q) (arr)[(size_t)(p) * (nso) + (q)]
#define IDX4(arr, nso, p, q, r, s)                                             \
  (arr)[(((size_t)(p) * (nso) + (q)) * (nso) + (r)) * (nso) + (s)]

static int spin_of(int p) { return p % 2; }
static int spatial_of(int p) { return p / 2; }

/*
 * Antisymmetrized physicist-notation spin-orbital integrals from the
 * spatial chemist-notation MO integrals:
 *   <pq|rs> = (pr|qs)                    [chemist -> physicist]
 *   <pq||rs> = <pq|rs> - <pq|sr>         [antisymmetrize]
 * with spin conservation (pq|rs)'s spin-orbital version is nonzero only when
 * spin_p=spin_r and spin_q=spin_s (for <pq|rs>), or spin_p=spin_s and
 * spin_q=spin_r (for <pq|sr>).
 */
static double *build_antisymmetrized(int n_spatial, const double *eri_mo,
                                     int nso) {
  double *v = calloc((size_t)nso * nso * nso * nso, sizeof(double));
  if (!v) {
    return NULL;
  }

  for (int p = 0; p < nso; p++) {
    int sp = spin_of(p), Pp = spatial_of(p);

    for (int q = 0; q < nso; q++) {
      int sq = spin_of(q), Pq = spatial_of(q);

      for (int r = 0; r < nso; r++) {
        int sr = spin_of(r), Pr = spatial_of(r);

        for (int s = 0; s < nso; s++) {
          int ss = spin_of(s), Ps = spatial_of(s);
          double v1 = (sp == sr && sq == ss)
                          ? MOLINT_ERI(eri_mo, n_spatial, Pp, Pr, Pq, Ps)
                          : 0.0;
          double v2 = (sp == ss && sq == sr)
                          ? MOLINT_ERI(eri_mo, n_spatial, Pp, Ps, Pq, Pr)
                          : 0.0;

          IDX4(v, nso, p, q, r, s) = v1 - v2;
        }
      }
    }
  }

  return v;
}

/* \tau_ij^ab = t2_ijab + t1_ia * t1_jb - t1_ib * t1_ja (full antisymmetrized
 * T1xT1 product, used in Wmnij/Wabef/T2-update). */
static double tau_full(const double *t1, const double *t2, int nso, int i,
                       int j, int a, int b) {
  return IDX4(t2, nso, i, j, a, b) + IDX2(t1, nso, i, a) * IDX2(t1, nso, j, b) -
         IDX2(t1, nso, i, b) * IDX2(t1, nso, j, a);
}

// \tau-\tilde_ij^ab = t2_ijab + 0.5*(t1_ia*t1_jb - t1_ib*t1_ja)
static double tau_tilde(const double *t1, const double *t2, int nso, int i,
                        int j, int a, int b) {
  return IDX4(t2, nso, i, j, a, b) +
         0.5 * (IDX2(t1, nso, i, a) * IDX2(t1, nso, j, b) -
                IDX2(t1, nso, i, b) * IDX2(t1, nso, j, a));
}

typedef struct {
  int nso;
  int nocc, nvirt;
  int *occ, *virt; /* nocc / nvirt long, spin-orbital indices */
  double *Fso;     /* nso x nso, diagonal (canonical orbitals) */
  double *V;       /* nso^4, antisymmetrized physicist integrals */
  double *Dia;     /* nso x nso, only occ x virt entries used */
  double *Dijab;   /* nso^4, only occ x occ x virt x virt entries used */
} ccsd_ctx_t;

static void build_intermediates(const ccsd_ctx_t *ctx, const double *t1,
                                const double *t2, double *Fae, double *Fmi,
                                double *Fme, double *Wmnij, double *Wabef,
                                double *Wmbej) {
  int nso = ctx->nso;
  const int *occ = ctx->occ, *virt = ctx->virt;
  int no = ctx->nocc, nv = ctx->nvirt;
  const double *Fso = ctx->Fso;
  const double *V = ctx->V;

  for (int ai = 0; ai < nv; ai++) {
    int a = virt[ai];

    for (int ei = 0; ei < nv; ei++) {
      int e = virt[ei];
      double s = 0.0;

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];
        s -= 0.5 * IDX2(Fso, nso, m, e) * IDX2(t1, nso, m, a);

        for (int fi = 0; fi < nv; fi++) {
          int f = virt[fi];
          s += IDX2(t1, nso, m, f) * IDX4(V, nso, m, a, f, e);
        }

        for (int ni = 0; ni < no; ni++) {
          int n = occ[ni];

          for (int fi = 0; fi < nv; fi++) {
            int f = virt[fi];
            s -= 0.5 * tau_tilde(t1, t2, nso, m, n, a, f) *
                 IDX4(V, nso, m, n, e, f);
          }
        }
      }

      IDX2(Fae, nso, a, e) = (a != e ? IDX2(Fso, nso, a, e) : 0.0) + s;
    }
  }

  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    for (int ii = 0; ii < no; ii++) {
      int i = occ[ii];
      double s = 0.0;

      for (int ei = 0; ei < nv; ei++) {
        int e = virt[ei];
        s += 0.5 * IDX2(t1, nso, i, e) * IDX2(Fso, nso, m, e);

        for (int ni = 0; ni < no; ni++) {
          int n = occ[ni];
          s += IDX2(t1, nso, n, e) * IDX4(V, nso, m, n, i, e);
        }

        for (int ni = 0; ni < no; ni++) {
          int n = occ[ni];

          for (int fi = 0; fi < nv; fi++) {
            int f = virt[fi];
            s += 0.5 * tau_tilde(t1, t2, nso, i, n, e, f) *
                 IDX4(V, nso, m, n, e, f);
          }
        }
      }

      IDX2(Fmi, nso, m, i) = (m != i ? IDX2(Fso, nso, m, i) : 0.0) + s;
    }
  }

  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    for (int ei = 0; ei < nv; ei++) {
      int e = virt[ei];
      double s = IDX2(Fso, nso, m, e);

      for (int ni = 0; ni < no; ni++) {
        int n = occ[ni];

        for (int fi = 0; fi < nv; fi++) {
          int f = virt[fi];
          s += IDX2(t1, nso, n, f) * IDX4(V, nso, m, n, e, f);
        }
      }

      IDX2(Fme, nso, m, e) = s;
    }
  }

  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    for (int ni = 0; ni < no; ni++) {
      int n = occ[ni];

      for (int ii = 0; ii < no; ii++) {
        int i = occ[ii];

        for (int ji = 0; ji < no; ji++) {
          int j = occ[ji];
          double s = IDX4(V, nso, m, n, i, j);

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];
            s += IDX2(t1, nso, j, e) * IDX4(V, nso, m, n, i, e) -
                 IDX2(t1, nso, i, e) * IDX4(V, nso, m, n, j, e);
          }

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];

            for (int fi = 0; fi < nv; fi++) {
              int f = virt[fi];
              s += 0.25 * tau_full(t1, t2, nso, i, j, e, f) *
                   IDX4(V, nso, m, n, e, f);
            }
          }

          IDX4(Wmnij, nso, m, n, i, j) = s;
        }
      }
    }
  }

  for (int ai = 0; ai < nv; ai++) {
    int a = virt[ai];

    for (int bi = 0; bi < nv; bi++) {
      int b = virt[bi];

      for (int ei = 0; ei < nv; ei++) {
        int e = virt[ei];

        for (int fi = 0; fi < nv; fi++) {
          int f = virt[fi];
          double s = IDX4(V, nso, a, b, e, f);

          for (int mi = 0; mi < no; mi++) {
            int m = occ[mi];
            s -= IDX2(t1, nso, m, b) * IDX4(V, nso, a, m, e, f) -
                 IDX2(t1, nso, m, a) * IDX4(V, nso, b, m, e, f);
          }

          for (int mi = 0; mi < no; mi++) {
            int m = occ[mi];

            for (int ni = 0; ni < no; ni++) {
              int n = occ[ni];
              s += 0.25 * tau_full(t1, t2, nso, m, n, a, b) *
                   IDX4(V, nso, m, n, e, f);
            }
          }

          IDX4(Wabef, nso, a, b, e, f) = s;
        }
      }
    }
  }

  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];

    for (int bi = 0; bi < nv; bi++) {
      int b = virt[bi];

      for (int ei = 0; ei < nv; ei++) {
        int e = virt[ei];

        for (int ji = 0; ji < no; ji++) {
          int j = occ[ji];
          double s = IDX4(V, nso, m, b, e, j);

          for (int fi = 0; fi < nv; fi++) {
            int f = virt[fi];
            s += IDX2(t1, nso, j, f) * IDX4(V, nso, m, b, e, f);
          }

          for (int ni = 0; ni < no; ni++) {
            int n = occ[ni];
            s -= IDX2(t1, nso, n, b) * IDX4(V, nso, m, n, e, j);
          }

          for (int ni = 0; ni < no; ni++) {
            int n = occ[ni];

            for (int fi = 0; fi < nv; fi++) {
              int f = virt[fi];
              s -= (0.5 * IDX4(t2, nso, j, n, f, b) +
                    IDX2(t1, nso, j, f) * IDX2(t1, nso, n, b)) *
                   IDX4(V, nso, m, n, e, f);
            }
          }

          IDX4(Wmbej, nso, m, b, e, j) = s;
        }
      }
    }
  }
}

static void update_amplitudes(const ccsd_ctx_t *ctx, const double *t1,
                              const double *t2, double *new_t1, double *new_t2,
                              double *Fae, double *Fmi, double *Fme,
                              double *Wmnij, double *Wabef, double *Wmbej) {
  int nso = ctx->nso;
  const int *occ = ctx->occ, *virt = ctx->virt;
  int no = ctx->nocc, nv = ctx->nvirt;
  const double *Fso = ctx->Fso;
  const double *V = ctx->V;

  build_intermediates(ctx, t1, t2, Fae, Fmi, Fme, Wmnij, Wabef, Wmbej);

  for (int ii = 0; ii < no; ii++) {
    int i = occ[ii];

    for (int ai = 0; ai < nv; ai++) {
      int a = virt[ai];
      double s = IDX2(Fso, nso, i, a);

      for (int ei = 0; ei < nv; ei++) {
        int e = virt[ei];
        s += IDX2(t1, nso, i, e) * IDX2(Fae, nso, a, e);
      }

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];
        s -= IDX2(t1, nso, m, a) * IDX2(Fmi, nso, m, i);
      }

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];

        for (int ei = 0; ei < nv; ei++) {
          int e = virt[ei];
          s += IDX4(t2, nso, i, m, a, e) * IDX2(Fme, nso, m, e);
        }
      }

      for (int ni = 0; ni < no; ni++) {
        int n = occ[ni];

        for (int fi = 0; fi < nv; fi++) {
          int f = virt[fi];
          s -= IDX2(t1, nso, n, f) * IDX4(V, nso, n, a, i, f);
        }
      }

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];

        for (int ei = 0; ei < nv; ei++) {
          int e = virt[ei];

          for (int fi = 0; fi < nv; fi++) {
            int f = virt[fi];
            s -= 0.5 * IDX4(t2, nso, i, m, e, f) * IDX4(V, nso, m, a, e, f);
          }
        }
      }

      for (int mi = 0; mi < no; mi++) {
        int m = occ[mi];

        for (int ei = 0; ei < nv; ei++) {
          int e = virt[ei];

          for (int ni = 0; ni < no; ni++) {
            int n = occ[ni];
            s -= 0.5 * IDX4(t2, nso, m, n, a, e) * IDX4(V, nso, n, m, e, i);
          }
        }
      }

      IDX2(new_t1, nso, i, a) = s / IDX2(ctx->Dia, nso, i, a);
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
          double s = IDX4(V, nso, i, j, a, b);

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];
            double corr_b = 0.0, corr_a = 0.0;

            for (int mi = 0; mi < no; mi++) {
              int m = occ[mi];

              corr_b += IDX2(t1, nso, m, b) * IDX2(Fme, nso, m, e);
              corr_a += IDX2(t1, nso, m, a) * IDX2(Fme, nso, m, e);
            }

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

            s -= IDX4(t2, nso, i, m, a, b) *
                     (IDX2(Fmi, nso, m, j) + 0.5 * corr_j) -
                 IDX4(t2, nso, j, m, a, b) *
                     (IDX2(Fmi, nso, m, i) + 0.5 * corr_i);
          }

          for (int mi = 0; mi < no; mi++) {
            int m = occ[mi];

            for (int ni = 0; ni < no; ni++) {
              int n = occ[ni];
              s += 0.5 * tau_full(t1, t2, nso, m, n, a, b) *
                   IDX4(Wmnij, nso, m, n, i, j);
            }
          }

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];

            for (int fi = 0; fi < nv; fi++) {
              int f = virt[fi];
              s += 0.5 * tau_full(t1, t2, nso, i, j, e, f) *
                   IDX4(Wabef, nso, a, b, e, f);
            }
          }

          for (int mi = 0; mi < no; mi++) {
            int m = occ[mi];

            for (int ei = 0; ei < nv; ei++) {
              int e = virt[ei];
              s += IDX4(t2, nso, i, m, a, e) * IDX4(Wmbej, nso, m, b, e, j) -
                   IDX2(t1, nso, i, e) * IDX2(t1, nso, m, a) *
                       IDX4(V, nso, m, b, e, j);
              s -= IDX4(t2, nso, i, m, b, e) * IDX4(Wmbej, nso, m, a, e, j) -
                   IDX2(t1, nso, i, e) * IDX2(t1, nso, m, b) *
                       IDX4(V, nso, m, a, e, j);
              s -= IDX4(t2, nso, j, m, a, e) * IDX4(Wmbej, nso, m, b, e, i) -
                   IDX2(t1, nso, j, e) * IDX2(t1, nso, m, a) *
                       IDX4(V, nso, m, b, e, i);
              s += IDX4(t2, nso, j, m, b, e) * IDX4(Wmbej, nso, m, a, e, i) -
                   IDX2(t1, nso, j, e) * IDX2(t1, nso, m, b) *
                       IDX4(V, nso, m, a, e, i);
            }
          }

          for (int ei = 0; ei < nv; ei++) {
            int e = virt[ei];
            s += IDX2(t1, nso, i, e) * IDX4(V, nso, a, b, e, j) -
                 IDX2(t1, nso, j, e) * IDX4(V, nso, a, b, e, i);
          }

          for (int mi = 0; mi < no; mi++) {
            int m = occ[mi];
            s -= IDX2(t1, nso, m, a) * IDX4(V, nso, m, b, i, j) -
                 IDX2(t1, nso, m, b) * IDX4(V, nso, m, a, i, j);
          }

          IDX4(new_t2, nso, i, j, a, b) = s / IDX4(ctx->Dijab, nso, i, j, a, b);
        }
      }
    }
  }
}

static double ccsd_energy(const ccsd_ctx_t *ctx, const double *t1,
                          const double *t2) {
  int nso = ctx->nso;
  const int *occ = ctx->occ, *virt = ctx->virt;
  int no = ctx->nocc, nv = ctx->nvirt;
  const double *Fso = ctx->Fso;
  const double *V = ctx->V;

  double E = 0.0;
  for (int ii = 0; ii < no; ii++) {
    int i = occ[ii];

    for (int ai = 0; ai < nv; ai++) {
      int a = virt[ai];

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
          double v = IDX4(V, nso, i, j, a, b);

          E += 0.25 * v * IDX4(t2, nso, i, j, a, b);
          E += 0.5 * v * IDX2(t1, nso, i, a) * IDX2(t1, nso, j, b);
        }
      }
    }
  }

  return E;
}

ccsd_result_t *ccsd_run(int n_spatial, const double *h_mo, const double *eri_mo,
                        const double *mo_energy, int n_electrons,
                        int n_frozen_spatial, double e_rhf, double conv_tol,
                        int max_iter) {
  if (n_spatial <= 0 || !h_mo || !eri_mo || !mo_energy ||
      n_electrons % 2 != 0 || 2 * n_frozen_spatial >= n_electrons) {
    return NULL;
  }

  (void)h_mo; /* canonical RHF: off-diagonal core Hamiltonian in the MO basis is
                 zero by Brillouin's theorem, only orbital energies (the
                 diagonal) are needed. */

  int nso = 2 * n_spatial;
  ccsd_ctx_t ctx;

  ctx.nso = nso;
  ctx.nocc = n_electrons - 2 * n_frozen_spatial;
  ctx.nvirt = nso - n_electrons;
  ctx.occ = malloc((size_t)ctx.nocc * sizeof(int));
  ctx.virt = malloc((size_t)ctx.nvirt * sizeof(int));
  ctx.Fso = calloc((size_t)nso * nso, sizeof(double));
  ctx.V = build_antisymmetrized(n_spatial, eri_mo, nso);
  ctx.Dia = calloc((size_t)nso * nso, sizeof(double));
  ctx.Dijab = calloc((size_t)nso * nso * nso * nso, sizeof(double));

  if (!ctx.occ || !ctx.virt || !ctx.Fso || !ctx.V || !ctx.Dia || !ctx.Dijab) {
    free(ctx.occ);
    free(ctx.virt);
    free(ctx.Fso);
    free(ctx.V);
    free(ctx.Dia);
    free(ctx.Dijab);

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

  for (int ii = 0; ii < no; ii++) {
    int i = ctx.occ[ii];

    for (int ji = 0; ji < no; ji++) {
      int j = ctx.occ[ji];

      for (int ai = 0; ai < nv; ai++) {
        int a = ctx.virt[ai];

        for (int bi = 0; bi < nv; bi++) {
          int b = ctx.virt[bi];

          IDX4(ctx.Dijab, nso, i, j, a, b) =
              IDX2(ctx.Fso, nso, i, i) + IDX2(ctx.Fso, nso, j, j) -
              IDX2(ctx.Fso, nso, a, a) - IDX2(ctx.Fso, nso, b, b);
        }
      }
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
  double *Wabef = calloc((size_t)nso * nso * nso * nso, sizeof(double));
  double *Wmbej = calloc((size_t)nso * nso * nso * nso, sizeof(double));

  ccsd_result_t *result = malloc(sizeof(ccsd_result_t));

  if (!t1 || !t2 || !new_t1 || !new_t2 || !Fae || !Fmi || !Fme || !Wmnij ||
      !Wabef || !Wmbej || !result) {
    free(t1);
    free(t2);
    free(new_t1);
    free(new_t2);
    free(Fae);
    free(Fmi);
    free(Fme);
    free(Wmnij);
    free(Wabef);
    free(Wmbej);
    free(result);
    free(ctx.occ);
    free(ctx.virt);
    free(ctx.Fso);
    free(ctx.V);
    free(ctx.Dia);
    free(ctx.Dijab);

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
              IDX4(ctx.V, nso, i, j, a, b) / IDX4(ctx.Dijab, nso, i, j, a, b);
        }
      }
    }
  }

  double e_prev = 0.0;
  int converged = 0, it;

  for (it = 0; it < max_iter; it++) {
    update_amplitudes(&ctx, t1, t2, new_t1, new_t2, Fae, Fmi, Fme, Wmnij, Wabef,
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

  free(t1);
  free(t2);
  free(new_t1);
  free(new_t2);
  free(Fae);
  free(Fmi);
  free(Fme);
  free(Wmnij);
  free(Wabef);
  free(Wmbej);
  free(ctx.occ);
  free(ctx.virt);
  free(ctx.Fso);
  free(ctx.V);
  free(ctx.Dia);
  free(ctx.Dijab);

  return result;
}
