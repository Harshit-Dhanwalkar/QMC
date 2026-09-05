#ifndef QMC_MP2_H
#define QMC_MP2_H

#include "hartree_fock.h"

/*
 * Second-order Moller-Plesset perturbation theory (MP2), restricted
 * closed-shell (Reference: Szabo & Ostlund Ch. 6).
 *
 * closed-shell RHF MP2, spatial orbitals, chemist notation:
 *   E_MP2 = \sum_{i,j occ} sum_{a,b virt} (ia|jb) *
 *           [2 * (ia|jb) - (ib|ja)] / (eps_i + eps_j - eps_a - eps_b)
 * Where
 *   (ia|jb) = \int \int \phi_i(1)\phi_a(1) (1/r12) \phi_j(2)\phi_b(2) dr1 dr2.
 *
 * Two implementations are provided:
 *  - mp2_correlation_energy(): restricted to s-only atomic orbitals
 *    (hartree_fock_atom_s_orbitals' radial-grid solver). For s-only orbitals,
 *    the multipole expansion of 1/r12 collapses to its L=0 term exactly, giving
 *    a closed form (ia|jb) = \int u_i(r) u_a(r) Y0_jb(r) dr (compute_Y0 kernel)
 *    - this only ever sees the l=0 slice  of the true MP2 correlation energy,
 *    since there are no p/d/... virtual orbitals to correlate into.
 *  - molecular_mp2(): general angular momentum, any basis set / molecule
 *    supported by molecular_hf.h's RHF solver (McMurchie-Davidson AO integrals
 *    + molecular_ao_to_mo, the same infrastructure ccsd_run() uses). This is
 *    what actually resolves the s-only restriction above - (ia|jb) here is a
 *    real p/d/...-including two-electron MO integral, not a multipole
 *    truncation.
 */

typedef struct {
  double e_hf;  // input HF total energy, Hartree
  double e_mp2; // MP2 correlation correction, Hartree (l=0-restricted, expected
                // -ve) */
  double e_total; // e_hf + e_mp2
  int n_occ;      // occupied orbitals used (= hf->n_orbitals)
  int n_virt;     // virtual orbitals used (<= hf->n_virtual)
} mp2_result_t;

/*
 * Compute (l=0-restricted) MP2 correlation energy on top of converged
 * hf_result_t.
 *
 * r, N     : radial grid passed to hartree_fock_atom_s_orbitals
 * n_virtual: number of lowest-energy virtual orbitals to include (1 <=
 *            n_virtual <= hf->n_virtual)
 *
 * Returns a zeroed mp2_result_t (e_mp2=0) on invalid input (NULL hf/r,
 * n_virtual out of range).
 */
mp2_result_t mp2_correlation_energy(const hf_result_t *hf, const double *r,
                                    int N, int n_virtual);

typedef struct {
  double e_rhf;   // input RHF total energy, Hartree
  double e_mp2;   // MP2 correlation correction, Hartree (expected <= 0)
  double e_total; // e_rhf + e_mp2
  int n_occ;      // occupied spatial orbitals actually correlated (excludes
                  // n_frozen_spatial)
  int n_virt;     // virtual spatial orbitals used (n_basis - n_electrons/2)
} molecular_mp2_result_t;

/*
 * General closed-shell RHF-MP2 correlation energy for any molecule/basis
 * set supported by molecular_hf.h - not restricted to s orbitals.
 *
 * n_basis, eri_mo, mo_energy: same convention as ccsd_run() - eri_mo is
 * n_basis^4 two-electron integrals in the MO basis, chemist notation
 * (MOLINT_ERI macro / molecular_ao_to_mo's output), mo_energy is n_basis RHF
 * orbital energies.
 *
 * Where
 *  - n_electrons     : total electron count (must be even: closed-shell)
 *  - n_frozen_spatial: number of lowest-energy occupied spatial orbitals to
 *                      exclude from the correlation sum (0 for none)
 *  - e_rhf           : converged RHF total energy, folded into e_total for
 *                      convenience
 *
 * Returns a zeroed molecular_mp2_result_t (e_mp2=0, n_occ=n_virt=0) on invalid
 * input: n_basis<=0, NULL eri_mo/mo_energy, n_electrons odd or out of [2 *
 * n_frozen_spatial + 2, 2 * n_basis] range.
 */
molecular_mp2_result_t molecular_mp2(int n_basis, const double *eri_mo,
                                     const double *mo_energy, int n_electrons,
                                     int n_frozen_spatial, double e_rhf);

#endif
