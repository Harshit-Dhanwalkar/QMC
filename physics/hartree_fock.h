#ifndef QMC_HARTREE_FOCK_H
#define QMC_HARTREE_FOCK_H

#include "../core/matrix.h"
#include "../core/vector.h"

/*
 * Restricted Hartree-Fock (RHF) self-consistent field for closed-shell
 * atoms/ions whose occupied subshells are all s-type (l = 0): He (1s^2),
 * Be (1s^2 2s^2), s-shell ions, etc. Atomic units throughout (\hbar = m_e =
 * e = 4 * \pi * eps0 = 1)
 *
 * NOTE: Only s-orbital: For central potential, e^--e^- interaction 1/|r-r'| is
 * expanded in Legendre multipoles Y^L. Angular (Gaunt-coefficient) selection
 * rules restrict which L actually contribute to direct (J) and exchange (K)
 * integrals between two orbitals of angular momenta l_a, l_b. For l_a = l_b = 0
 * only surviving multipole is L = 0, so both J and K collapse onto same radial
 * kernel 1/max(r,r'), and exchange operator becomes an exact, closed radial
 * expression whihc discretize directly.
 *
 * MODEL:
 *   n_orbitals distinct radial s-orbitals u_1(r) < u_2(r) < ... (u = r*R,
 *   normalized: integral u_k(r)^2 dr = 1), each 'doubly occupied' (Aufbau, spin
 *   up + spin down), i.e. this models atoms/ions with e^- count = 2 *
 *   n_orbitals and all e^-s in s subshells (He: n_orbitals=1, Be: n_orbitals=2,
 *   ...).
 *
 * Closed-shell RHF Fock operator (Refrence: Szabo & Ostlund eq. 3.184):
 *   F = h + \sum_k ( 2 * J_k - K_k ),
 *   h = -(1/2) * d^2/dr^2 + l(l + 1) / (2 * r^2) - Z/r
 * discretized as a dense real-symmetric matrix, diagonalized each SCF iteration
 * via Hermitian eigensolver.
 *
 * Total electronic energy (Szabo & Ostlund eq. 3.184):
 *   E = \sum_k 2 * eps_k - \sum_{i,j} (2 J_ij - K_ij)
 * Where
 *  eps_k are converged orbital (Fock) eigenvalues and J_ij, K_ij are
 * direct/exchange two-e^-s integrals between orbitals i, j.
 */

typedef struct {
  int n_orbitals;           // number of doubly-occupied s-orbitals
  int N;                    // radial grid size (same N as input r[])
  double Z;                 // nuclear charge used
  double *orbital_energies; // size n_orbitals, converged Fock eigenvalues
  cvector_t **orbitals;     // size n_orbitals; orbitals[k] has length N,
                            // u_k(r) = r * R_k(r) stored in .re part,
                            // normalized: \int u_k(r)^2 dr = 1
  double total_energy;      // converged total electronic energy (Hartree)
  int iterations;           // SCF iterations actually performed
  int converged;            // 1 if converged within max_iter, else 0

  /* NOTE: Unoccupied ("virtual") Fock eigenpairs from converged Fock matrix
   * diagonalization n_virtual = N - n_orbitals. convention/normalization as
   * occupied `orbitals` above, ascending energy order. */
  int n_virtual;
  double *virtual_energies;     // size n_virtual
  cvector_t **virtual_orbitals; // size n_virtual
} hf_result_t;

/*
 * Run RHF SCF for a closed-shell, s-orbitals-only atom/ion.
 *
 * r: uniform radial grid r[0..N-1], r[0] > 0 recommended (avoids r=0
 *    Coulomb singularity), r[1]-r[0] = dr constant (same convention as
 *    central_potential_radial_solve).
 * N: grid size, N >= 10 or so for anything meaningful.
 * Z: nuclear charge (e.g. 2.0 for helium, 4.0 for beryllium).
 * n_orbitals: number of distinct doubly-occupied s subshells to fill
 *    (1 = He-like/1s^2, 2 = Be-like/1s^2 2s^2, ...).
 * mix: linear density-mixing fraction in (0,1] applied to each SCF update
 *    (new = (1-mix)*old + mix*new_raw) for stability; 1.0 = no damping,
 *    smaller values (e.g. 0.3-0.5) trade convergence speed for stability.
 * tol: convergence threshold on max pointwise orbital change b/w successive
 *    (mixed) iterations.
 * max_iter: maximum SCF iterations before giving up (result->converged=0).
 *
 * Returns a allocated hf_result_t, or NULL on invalid input / allocation
 * failure.
 */
hf_result_t *hartree_fock_atom_s_orbitals(const double *r, int N, double Z,
                                          int n_orbitals, double mix,
                                          double tol, int max_iter);

void hf_result_free(hf_result_t *res);

/*
 * l=0 ("monopole") radial Coulomb kernel shared with post-HF methods: L=0 term
 * of 1/r12 multipole expansion, i.e. electrostatic potential at r of a
 * spherical charge density proportional to u_a(r')u_b(r'):
 *   Y0_ab(r) = (1/r) * \int_0^r u_a(r')u_b(r') dr' +
 *              \int_r^r_{max} [u_a(r')u_b(r') / r'] dr'
 * Writes N values into Y0_out (length N).
 */
void compute_Y0(const double *r, int N, double dr, const double *ua,
                const double *ub, double *Y0_out);

#endif // QMC_HARTREE_FOCK_H
