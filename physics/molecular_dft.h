#ifndef QMC_MOLECULAR_DFT_H
#define QMC_MOLECULAR_DFT_H

/*
 * Molecular (multi-atom) Kohn-Sham LDA DFT via a Becke fuzzy-Voronoi numerical
 * integration grid.
 *
 * NOTE: dft.c already implements atomic radial-grid KS-LDA (single center, uses
 * spherical symmetry to reduce problem to 1D). This module generalizes to
 * arbitrary molecules with no special symmetry, reusing dft.c's pointwise LDA
 * exchange-correlation functions ( Reference: lda_xc_energy_density,
 * lda_xc_potential : Slater exchange + Perdew-Zunger 1981 correlation) but
 * evaluated over a 3D grid built from Becke's atom-centered fuzzy-cell
 * partition (References: Becke 1988, J. Chem. Phys. 88, 2547):
 *
 *   1. Each atom gets its own atom-centered (radial x angular) quadrature
 *      grid covering all of space.
 *   2. Becke's iterated-polynomial "fuzzy Voronoi" weight function w_A(r)
 *      smoothly partitions space between atoms (w_A(r) -> 1 near atom A's
 *      nucleus, -> 0 near any other atom, with a smooth transition region
 *      in between) so that summing every atom's weighted grid gives a
 *      partition of unity: sum_A w_A(r) = 1 everywhere.
 *   3. The molecular integral of any function f(r) is then
 *      sum_A sum_{g in atom A's grid} w_A(r_g) * quad_weight_g * f(r_g).
 *
 * Where
 *  - Radial grid: Becke's Gauss-Chebyshev-of-the-second-kind mapping, r_i =
 *                 r_m*(1 * + x_i) / (1 - x_i), x_i = \cos(i * \pi / (N + 1)),
 *                 which maps (-1,1) to (0, \infty).
 *  - Angular grid: a Gauss-Legendre (in \cos(\theta)) x uniform-trapezoidal
 *                 (in \phi) product grid but exact for its quadrature order.
 */

#include "molecular_integrals.h"

typedef struct {
  double x, y, z;
  double weight; /* combined radial * angular * Jacobian * Becke-atomic-
                  * partition-weight quadrature weight at this point */
} dft_grid_point_t;

typedef struct {
  dft_grid_point_t *points;
  int n_points;
} molecular_grid_t;

/*
 * Build the full molecular Becke grid: n_radial points per atom's radial
 * shell, n_polar x n_azimuthal points per angular shell (so each atom
 * contributes n_radial*n_polar*n_azimuthal points, before any are effectively
 * zeroed out by a near-zero Becke weight far from that atom). radial_scale sets
 * Gauss-Chebyshev radial mapping's r_m parameter (in bohr) : general-purpose
 * default is 1.0.
 */
molecular_grid_t *molecular_grid_build(const molecule_t *mol, int n_radial,
                                       int n_polar, int n_azimuthal,
                                       double radial_scale);

/* Wrapper: molecular_grid_build with parameters
 * HACK: Larger/more diffuse molecules may need a denser grid; pass explicit
 * parameters to molecular_grid_build directly for that case. */
molecular_grid_t *molecular_grid_build_default(const molecule_t *mol);

void molecular_grid_free(molecular_grid_t *grid);

typedef struct {
  double total_energy;
  double e_core;            /* Tr(D*Hcore) */
  double e_coulomb;         /* (1/2) * Tr(D * J), the classical (Hartree)
                               electron-electron repulsion energy */
  double e_xc;              /* integral n(r) * eps_xc(n(r)) dr */
  double e_nuclear;         /* classical nuclear-nuclear repulsion */
  double *orbital_energies; /* length n_basis, ascending */
  cmatrix_t *C;             /* MO coefficients, n_basis x n_basis */
  int n_basis;
  int n_electrons;
  int converged;
  int iterations;
} molecular_dft_result_t;

/*
 * Restricted (closed-shell) Kohn-Sham LDA SCF for a general molecule/basis.
 * Reuses molecular_integrals.c for S/Hcore/ERI (so classical Coulomb (J) term
 * is exact, analytic : only exchange-correlation term is evaluated numerically
 * on `grid`) and dft.c's lda_xc_energy_density / lda_xc_potential for the
 * pointwise LDA functional.
 *
 * NOTE: `mix` is a linear density-mixing/damping factor in (0,1]:
 *   D_new = mix * D_computed + (1 - mix) * D_old.
 * Plain unmixed SCF (mix=1.0) converges for simplest systems (e.g. H2) but
 * oscillates for anything with near-degenerate orbitals (e.g. LiH), mix=0.3 is
 * a default
 */
molecular_dft_result_t *molecular_ks_lda(basis_function_t **basis, int n_basis,
                                         const molecule_t *mol, int n_electrons,
                                         const molecular_grid_t *grid,
                                         double mix, double conv_tol,
                                         int max_iter);

/* molecular_ks_lda with mix=0.3, conv_tol=1e-9, max_iter=200 */
molecular_dft_result_t *molecular_ks_lda_default(basis_function_t **basis,
                                                 int n_basis,
                                                 const molecule_t *mol,
                                                 int n_electrons,
                                                 const molecular_grid_t *grid);

/*
 * Restricted (closed-shell) Kohn-Sham PBE (GGA) SCF: same interface and
 * approach as molecular_ks_lda, but uses dft.h's pbe_xc_energy_density /
 * pbe_xc_potential in place of LDA functional, which additionally needs the
 * density gradient (not just density) at every grid point. That requires each
 * basis function's gradient (molecular_integrals.h's basis_function_gradient),
 * not just its value, and GGA Kohn-Sham potential matrix element
 *   V_pq = \int [ vrho * \phi_p * \phi_q +
 *                 2 * vsigma*(\grad n).(\phi_q * \grad(\phi_p) +
 *                 \phi_p * \grad(\phi_q)) ] dr
 *
 * Where
 *   vrho = d(n*eps_xc) / dn
 *   vsigma = d(n * eps_xc) / dsigma,
 *   \sigma=|grad n|^2),
 * in place of LDA's purely-diagonal-in-density V_pq = \int[vxc * \phi_p *
 * \phi_q]dr.
 */
molecular_dft_result_t *molecular_ks_pbe(basis_function_t **basis, int n_basis,
                                         const molecule_t *mol, int n_electrons,
                                         const molecular_grid_t *grid,
                                         double mix, double conv_tol,
                                         int max_iter);

/* Convenience wrapper: molecular_ks_pbe with mix=0.3, conv_tol=1e-9,
 * max_iter=200 */
molecular_dft_result_t *molecular_ks_pbe_default(basis_function_t **basis,
                                                 int n_basis,
                                                 const molecule_t *mol,
                                                 int n_electrons,
                                                 const molecular_grid_t *grid);

void molecular_dft_result_free(molecular_dft_result_t *res);

#endif
