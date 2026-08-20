#include "ccsd_t.h"
#include <stdlib.h>

#define IDX2(arr, nso, p, q) (arr)[(size_t)(p) * (nso) + (q)]
#define IDX4(arr, nso, p, q, r, s)                                             \
  (arr)[(((size_t)(p) * (nso) + (q)) * (nso) + (r)) * (nso) + (s)]

/* tmp_d(i,j,k,a,b,c) = t1_ia * <jk||bc> -- the disconnected-triples "base"
   term (before the P(i/jk)P(a/bc) antisymmetrization below is applied). */
static double tmp_d_val(const ccsd_amplitudes_t *amp, int i, int j, int k,
                        int a, int b, int c) {
  (void)k;
  int nso = amp->nso;
  return IDX2(amp->t1, nso, i, a) * IDX4(amp->V, nso, j, k, b, c);
}

/* tmp_c(i,j,k,a,b,c) = sum_e t2_jkae*<ei||bc> - sum_m t2_imbc*<majk> -- the
   connected-triples "base" term. */
static double tmp_c_val(const ccsd_amplitudes_t *amp, int i, int j, int k,
                        int a, int b, int c) {
  int nso = amp->nso;
  int nv = amp->nvirt, no = amp->nocc;
  const int *virt = amp->virt;
  const int *occ = amp->occ;

  double s = 0.0;
  for (int ei = 0; ei < nv; ei++) {
    int e = virt[ei];
    s += IDX4(amp->t2, nso, j, k, a, e) * IDX4(amp->V, nso, e, i, b, c);
  }
  for (int mi = 0; mi < no; mi++) {
    int m = occ[mi];
    s -= IDX4(amp->t2, nso, i, m, b, c) * IDX4(amp->V, nso, m, a, j, k);
  }
  return s;
}

typedef double (*tmp_fn_t)(const ccsd_amplitudes_t *, int, int, int, int, int,
                           int);

/* P(i/jk) P(a/bc) tmp(i,j,k,a,b,c): the 3-term cyclic antisymmetrizer over
   {i swapped with j, i swapped with k} composed with the same over
   {a swapped with b, a swapped with c} -- 9 terms total (3x3), signs from
   expanding (1 - Sij - Sik)(1 - Sab - Sac). */
static double P_ijk_abc(tmp_fn_t tmp, const ccsd_amplitudes_t *amp, int i,
                        int j, int k, int a, int b, int c) {
  double t = 0.0;
  t += tmp(amp, i, j, k, a, b, c);
  t -= tmp(amp, i, j, k, b, a, c);
  t -= tmp(amp, i, j, k, c, b, a);
  t -= tmp(amp, j, i, k, a, b, c);
  t += tmp(amp, j, i, k, b, a, c);
  t += tmp(amp, j, i, k, c, b, a);
  t -= tmp(amp, k, j, i, a, b, c);
  t += tmp(amp, k, j, i, b, a, c);
  t += tmp(amp, k, j, i, c, b, a);
  return t;
}

static double perturbative_triples(const ccsd_amplitudes_t *amp) {
  int no = amp->nocc, nv = amp->nvirt, nso = amp->nso;
  const int *occ = amp->occ, *virt = amp->virt;
  const double *Fso = amp->Fso;

  double E = 0.0;

  for (int ii = 0; ii < no; ii++) {
    int i = occ[ii];
    for (int ji = 0; ji < no; ji++) {
      int j = occ[ji];
      for (int ki = 0; ki < no; ki++) {
        int k = occ[ki];

        double fi = IDX2(Fso, nso, i, i);
        double fj = IDX2(Fso, nso, j, j);
        double fk = IDX2(Fso, nso, k, k);

        for (int ai = 0; ai < nv; ai++) {
          int a = virt[ai];
          for (int bi = 0; bi < nv; bi++) {
            int b = virt[bi];
            for (int ci = 0; ci < nv; ci++) {
              int c = virt[ci];

              double fa = IDX2(Fso, nso, a, a);
              double fb = IDX2(Fso, nso, b, b);
              double fc = IDX2(Fso, nso, c, c);
              double D = fi + fj + fk - fa - fb - fc;

              double t3d = P_ijk_abc(tmp_d_val, amp, i, j, k, a, b, c);
              double t3c = P_ijk_abc(tmp_c_val, amp, i, j, k, a, b, c);

              E += t3c * (t3c + t3d) / D;
            }
          }
        }
      }
    }
  }

  return E / 36.0;
}

ccsdt_result_t *ccsdt_run(int n_spatial, const double *h_mo,
                          const double *eri_mo, const double *mo_energy,
                          int n_electrons, int n_frozen_spatial, double e_rhf,
                          double conv_tol, int max_iter) {
  ccsd_amplitudes_t *amp = NULL;
  ccsd_result_t *ccsd =
      ccsd_run_ex(n_spatial, h_mo, eri_mo, mo_energy, n_electrons,
                  n_frozen_spatial, e_rhf, conv_tol, max_iter, &amp);

  if (!ccsd) {
    return NULL;
  }

  if (!ccsd->converged) {
    free(ccsd);
    ccsd_amplitudes_free(amp);
    return NULL;
  }

  double pert_t = perturbative_triples(amp);

  ccsdt_result_t *result = malloc(sizeof(ccsdt_result_t));
  if (result) {
    result->ccsd_correlation_energy = ccsd->correlation_energy;
    result->perturbative_correction = pert_t;
    result->total_energy = ccsd->total_energy + pert_t;
    result->ccsd_converged = ccsd->converged;
    result->ccsd_iterations = ccsd->iterations;
  }

  free(ccsd);
  ccsd_amplitudes_free(amp);

  return result;
}
