#ifndef QMC_CCSD_H
#define QMC_CCSD_H

/*
 * Closed-shell-reference spin-orbital CCSD (Coupled Cluster Singles and
 * Doubles), built on molecular_hf.c's canonical RHF orbitals and
 * molecular_ao_to_mo's MO-basis integrals.
 *
 *  NOTE: Uses standard spin-orbital formulation (Reference: Stanton, Gauss,
 * Watts & Bartlett, J. Chem. Phys. 94, 4334 (1991)): antisymmetrized
 * spin-orbital integrals built from the spatial RHF MOs (doubled into
 * alternating alpha/beta spin orbitals), then the T1/T2 amplitude equations
 * solved by plain fixed-point iteration through Fae/Fmi/Fme/Wmnij/Wabef/Wmbej
 * intermediates.
 *  HACK: Spin-orbital (not spin-adapted spatial-orbital) CCSD is used because
 * its equations are simpler to state and verify correctly.
 *  HACK: acceptable for
 * currently small test molecules (H2, LiH)
 *  TODO: scale to large systems.
 */

#include "molecular_hf.h"

typedef struct {
  double correlation_energy; /* E_CCSD - E_RHF */
  double total_energy;       /* E_RHF + correlation_energy */
  int converged;
  int iterations;
} ccsd_result_t;

/*
 * Converged CCSD amplitudes and underlying spin-orbital machinery they
 * were computed with, retained (rather than freed internally, as ccsd_run
 * does) for downstream use.
 * The perturbative (T) triples correction (ccsd_t.c), which needs T1/T2 plus
 * the antisymmetrized spin-orbital integrals V and Fock diagonal Fso to build
 * the connected / disconnected T3 contractions. Same flat-array indexing
 * convention as ccsd.c internals: IDX2(arr, nso, p, q) = arr[p * nso + q], IDX4
 * similarly; occ[]/virt[] list the nocc/nvirt active (frozen-core-excluded)
 * spin-orbital indices actually populated in Fso/V/t1/t2 (all are allocated at
 * full nso/nso^2/nso^4 size, but only entries touching occ[]/virt[] indices are
 * meaningful).
 */
typedef struct {
  int nso, nocc, nvirt;
  int *occ, *virt;
  double *Fso;
  double *V;
  double *t1;
  double *t2;
} ccsd_amplitudes_t;

void ccsd_amplitudes_free(ccsd_amplitudes_t *amp);

/*
 * Runs CCSD on top of a converged RHF result.
 *
 * n_spatial  : number of spatial MOs (basis functions).
 * h_mo       : n_spatial x n_spatial core Hamiltonian in the MO basis (from
 *              molecular_ao_to_mo).
 * eri_mo     : n_spatial^4 flat MO-basis ERI tensor, chemist notation (pq|rs),
 *              indexed via MOLINT_ERI (from molecular_ao_to_mo).
 * mo_energy  : n_spatial canonical RHF orbital energies (e.g.
 *              molecular_hf_result_t::orbital_energies).
 * n_electrons: total electron count (must be even; closed-shell reference).
 * n_frozen_spatial: number of lowest-energy spatial orbitals to exclude from
 *              the correlation treatment (0 for none). Must satisfy
 *              2*n_frozen_spatial < n_electrons. e_rhf: converged RHF total
 *              energy (electronic + nuclear repulsion), e.g.
 *              molecular_hf_result_t::total_energy.
 * conv_tol   : convergence threshold on the correlation energy between
 *              iterations.
 * max_iter   : iteration cap.
 *
 * Returns a newly allocated ccsd_result_t, or NULL on invalid input (odd
 * electron count, n_frozen_spatial too large, non-positive n_spatial) or
 * allocation failure.
 */
ccsd_result_t *ccsd_run(int n_spatial, const double *h_mo, const double *eri_mo,
                        const double *mo_energy, int n_electrons,
                        int n_frozen_spatial, double e_rhf, double conv_tol,
                        int max_iter);

/*
 * Same as ccsd_run, but if `amplitudes_out` is non-NULL, also returns converged
 * T1/T2 amplitudes and supporting spin-orbital machinery instead of freeing
 * them internally. Passing NULL for `amplitudes_out` makes this behave
 * identically to ccsd_run.
 */
ccsd_result_t *ccsd_run_ex(int n_spatial, const double *h_mo,
                           const double *eri_mo, const double *mo_energy,
                           int n_electrons, int n_frozen_spatial, double e_rhf,
                           double conv_tol, int max_iter,
                           ccsd_amplitudes_t **amplitudes_out);

#endif
