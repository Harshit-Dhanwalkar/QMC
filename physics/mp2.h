#ifndef QMC_MP2_H
#define QMC_MP2_H

#include "hartree_fock.h"

/*
 * Second-order Moller-Plesset perturbation theory (MP2), restricted
 * closed-shell (Refermce : Szabo & Ostlund Ch. 6).
 *
 * NOTE: : Since underlying HF is s-orbitals-only (occupied and virtual orbitals
 * here are all l=0), this computes only l=0-restricted slice of true MP2
 * correlation energy
 * TODO: Include double excitations into p, d, ... virtual orbitals (angular
 * correlation)
 *
 * closed-shell RHF MP2, spatial orbitals, chemist notation:
 *   E_MP2 = \sum_{i,j occ} sum_{a,b virt} (ia|jb) *
 *           [2 * (ia|jb) - (ib|ja)] / (eps_i + eps_j - eps_a - eps_b)
 * Where
 *   (ia|jb) = \int \int \phi_i(1)\phi_a(1) (1/r12) \phi_j(2)\phi_b(2) dr1 dr2.
 *
 * For s-only orbitals, multipole expansion of 1/r12 collapses to its L=0 term
 * exactly, giving closed form (ia|jb) = \int u_i(r) u_a(r) Y0_jb(r) dr
 * compute_Y0 kernel.
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
 * r, N     : radial grid passed to hartree_fock_atom_s_orbitals.
 * n_virtual: number of lowest-energy virtual orbitals to include (1 <=
 *            n_virtual <= hf->n_virtual).
 *
 * Returns a zeroed mp2_result_t (e_mp2=0) on invalid input (NULL hf/r,
 * n_virtual out of range).
 */
mp2_result_t mp2_correlation_energy(const hf_result_t *hf, const double *r,
                                    int N, int n_virtual);

#endif
