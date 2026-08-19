#ifndef QMC_CCSD_H
#define QMC_CCSD_H

/*
 * Closed-shell-reference spin-orbital CCSD (Coupled Cluster Singles and
 * Doubles), built on molecular_hf.c's canonical RHF orbitals and
 * molecular_ao_to_mo's MO-basis integrals.
 *
 * NOTE: Uses standard spin-orbital formulation (Reference: Stanton, Gauss,
 * Watts & Bartlett, J. Chem. Phys. 94, 4334 (1991)): antisymmetrized
 * spin-orbital integrals built from the spatial RHF MOs (doubled into
 * alternating alpha/beta spin orbitals), then the T1/T2 amplitude equations
 * solved by plain fixed-point iteration through the
 * Fae/Fmi/Fme/Wmnij/Wabef/Wmbej intermediates. Spin-orbital (not spin-adapted
 * spatial-orbital) CCSD is used because its equations are simpler to state and
 * verify correctly (fewer distinct spin cases to track).
 * // HACK: acceptable for currently small test molecules (H2, LiH)
 * // TODO: scale to large systems.
 */

#include "molecular_hf.h"

typedef struct {
  double correlation_energy; /* E_CCSD - E_RHF */
  double total_energy;       /* E_RHF + correlation_energy */
  int converged;
  int iterations;
} ccsd_result_t;

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

#endif
