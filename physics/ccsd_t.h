#ifndef QMC_CCSD_T_H
#define QMC_CCSD_T_H

/*
 * Perturbative triples correction, CCSD(T) (Reference: Raghavachari, Trucks,
 * Pople & Head-Gordon, Chem. Phys. Lett. 157, 479 (1989)), converged T1/T2
 * amplitudes (via ccsd_run_ex).
 *
 * Spin-orbital formulas (Crawford & Schaefer's "An Introduction to Coupled
 * Cluster Theory", project #6):
 *
 *   t3d_ijkabc = P(i / jk) P(a / bc) [ t1_ia * <jk||bc> ]
 *   t3c_ijkabc = P(i / jk) P(a / bc) [ \sum_e t2_jkae * <ei||bc>
 *                                    - \sum_m t2_imbc * <majk> ]
 *   E_(T) = (1/36) * \sum_ijkabc t3c_ijkabc * (t3c_ijkabc + t3d_ijkabc)
 *                                              / D_ijkabc
 *
 * Where :
 *  P(p/qr) f(p,q,r) = f(p,q,r) - f(q,p,r) - f(r,q,p)
 * the 3-term cyclic antisymmetrizer over the index being singled out), and
 * D_ijkabc is the usual triples energy denominator
 * f_ii + f_jj + f_kk - f_aa - f_bb - f_cc.
 */

#include "ccsd.h"

typedef struct {
  double ccsd_correlation_energy; /* E_CCSD - E_RHF, (same as ccsd_result_t) */
  double perturbative_correction; /* (T) term */
  double total_energy;            /* E_RHF + CCSD corr + (T) */
  int ccsd_converged;
  int ccsd_iterations;
} ccsdt_result_t;

/*
 * Same inputs as ccsd_run : runs CCSD to convergence, then adds the
 * perturbative (T) correction using the converged amplitudes.
 *
 * Returns NULL on the same invalid-input conditions as ccsd_run, or if the
 * underlying CCSD does not converge within max_iter
 *
 * HACK: ccsdt_result_t's ccsd_converged field lets the caller distinguish
 * "didn't converge" from "invalid input" if a non-NULL-on-non-convergence
 * variant is ever needed : for now this returns NULL in both cases, matching
 * ccsd_run's contract of not returning results from an unconverged
 * calculation.
 */
ccsdt_result_t *ccsdt_run(int n_spatial, const double *h_mo,
                          const double *eri_mo, const double *mo_energy,
                          int n_electrons, int n_frozen_spatial, double e_rhf,
                          double conv_tol, int max_iter);

#endif
