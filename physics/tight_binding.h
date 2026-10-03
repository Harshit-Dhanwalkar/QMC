#ifndef QMC_TIGHT_BINDING_H
#define QMC_TIGHT_BINDING_H

#include "../core/complex.h"
#include "../core/matrix.h"
#include <limits.h>

/*
 * Generic multi-band tight-binding models on 1D and 2D Bravais lattices:
 * Bloch Hamiltonians, band structures, density of states, and topological
 * invariants of occupied bands (Berry curvature, Chern number, Wilson-loop
 * / Zak phase)
 *
 * Natural units (\hbar = 1). Wave vectors are in REDUCED coordinates,
 * (k1, k2) = (k . a1, k . a2), so Brillouin zone is [0, 2 \pi)^2 for every
 * lattice and no geometry is needed. Physical velocities and effective masses
 * require lattice vectors and are left to caller
 *
 * NOTE: Model. A unit cell holds n_orb orbitals, all sitting at cell
 * origin ("periodic gauge"). Each hopping term is
 *   H += t * c^\dagger_{i,R} c_{j,R + n1 a1 + n2 a2} + h.c.
 * (added with tb_add_hopping), and each orbital has an on-site energy
 * The Bloch Hamiltonian is then
 *   H_ij(k) += t \xpe^{i (k1 n1 + k2 n2)}
 *   H_ji(k) += conj(t e^{i (k1 n1 + k2 n2)})
 *   H_ii(k)  = eps_i + ... and is exactly periodic, H(k + 2 \pi) = H(k),
 * which is what invariants below need. Because orbitals are placed at
 * cell origin, Berry phases are those of "periodic gauge" (Zak phase is
 * then quantized to 0 or \pi for inversion-symmetric chains, with no dependence
 * on where cell is cut)
 *
 * NOTE: Invariants use gauge-invariant discretized (Reference:
 * Fukui-Hatsugai-Suzuki 2005) method: eigenvectors on an n_k-point uniform grid
 * are compared only through determinants of overlap matrices det
 * <u_a(k)|u_b(k')>, so no smooth gauge is needed, band groups (several bands
 * treated together, non-Abelian Berry phase) are handled with same code, and a
 * group that is gapped from rest of spectrum gives an integer-quantized Chern
 * number as soon as grid resolves Berry curvature Berry flux of a plaquette:
 * - F = arg(U1 U2 U3 U4), product of four link variables U =
 *       det<u|u'>/|det<u|u'>| around plaquette traversed counter-clockwise in
 *       (k1, k2); each F lies in (-\pi, \pi]
 * - Chern number: C = (1/2 \pi) \sum F
 * - Zak / Wilson phase: \gamma = -arg \prod_k det<u_k|u_{k+1}> along k1
 * - Bands are numbered from 0 in ascending energy order Sign of C is a
 *        convention (orientation of k1 -> k2 plane and sign of Berry
 *        connection)
 *
 *  This module's convention is fixed by analytic references below and is shared
 * with lattice_hofstadter_chern_numbers() in lattice.h.
 */

typedef struct tb_model tb_model_t;

typedef enum {
  TB_OK = 0,
  TB_ERR_INVALID = -1,   // bad argument
  TB_ERR_MEMORY = -2,    // allocation failure
  TB_ERR_DEGENERATE = -3 // selected bands touch their neighbors, or the grid
                         // is too coarse (a link determinant vanishes)
} tb_status_t;

#define TB_MAX_DIM 2
#define TB_MAX_ORBITALS 64
#define TB_MAX_GRID 2048
#define TB_MAX_HOP_RANGE 64
#define TB_CHERN_UNDEFINED INT_MIN

/*
 * Create an empty model: `dim` = 1 or 2 spatial dimensions and `n_orb` orbitals
 * per cell (1 .. TB_MAX_ORBITALS), all on-site energies zero, no hoppings
 *
 * Returns NULL on invalid input or allocation failure
 */
tb_model_t *tb_model_alloc(int dim, int n_orb);
void tb_model_free(tb_model_t *model);

// Spatial dimension / orbital count, or 0 for a NULL model
int tb_model_dim(const tb_model_t *model);
int tb_model_orbitals(const tb_model_t *model);

tb_status_t tb_set_onsite(tb_model_t *model, int orbital, double energy);

/*
 * Add hopping t c^\dagger_{i,R} c_{j,R + n1 a1 + n2 a2} + h.c.
 *
 * i, j       : orbital indices
 * n1, n2     : cell offset (|n| <= TB_MAX_HOP_RANGE; n2 must be 0 when
 *              dim == 1)
 * t          : complex amplitude (finite)
 * An on-cell hopping of an orbital to itself (i == j, n1 == n2 == 0) is an
 * on-site energy and is rejected: use tb_set_onsite. Hoppings accumulate, so
 * adding same bond twice doubles it. A hopping i == j with (n1, n2) != 0 is
 * allowed and contributes 2 Re(t \exp^{i k.n}) to H_ii
 */
tb_status_t tb_add_hopping(tb_model_t *model, int i, int j, int n1, int n2,
                           complex_t t);

/*
 * Bloch Hamiltonian H(k1, k2) as a newly allocated n_orb x n_orb Hermitian
 * matrix (caller frees with cmatrix_free), or NULL on invalid input
 * For a 1D model k2 is ignored
 */
cmatrix_t *tb_hamiltonian(const tb_model_t *model, double k1, double k2);

/*
 * The n_orb band energies at k, ascending, written to `energies` (caller-
 * allocated, length n_orb)
 */
tb_status_t tb_bands(const tb_model_t *model, double k1, double k2,
                     double *energies);

/*
 * Band energies at n_points (>= 2) points evenly spaced on straight line
 * from (k1a, k2a) to (k1b, k2b), endpoints included. `energies_out` has length
 * n_points * n_orb and is laid out as energies_out[point * n_orb + band]
 */
tb_status_t tb_bands_along_line(const tb_model_t *model, double k1a, double k2a,
                                double k1b, double k2b, int n_points,
                                double *energies_out);

/*
 * Density of states with Gaussian broadening of width sigma, from a uniform
 * k-grid (n_k points for 1D, n_k x n_k for 2D, both starting at k = 0):
 *   dos(E) = (1 / N_k) \sum_{k, n} exp(-(E - E_n(k))^2 / 2 \sigma^2) /
 *            (\sigma \sqrt{2 \pi})
 * integrates to n_orb over whole band range (states per unit cell)
 * Evaluated on n_e (>= 2) energies evenly spaced from e_min to e_max
 * (inclusive); `dos_out` has length n_e
 */
tb_status_t tb_dos(const tb_model_t *model, int n_k, double e_min, double e_max,
                   int n_e, double sigma, double *dos_out);

/*
 * Berry flux through every plaquette of n_k x n_k grid for band group
 * [first_band, first_band + n_bands)
 * 2D models only
 *
 * flux_out : caller-allocated, length n_k * n_k; flux_out[a * n_k + b] is
 *            flux of plaquette whose lower-left corner is
 *            k = (2 \pi a / n_k, 2 \pi b / n_k). Its \sum / 2 \pi is Chern
 *            number; its pattern shows where Berry curvature lives (e.g. Dirac
 *            points of a gapped graphene-like model)
 * min_gap  : optional (may be NULL); receives smallest gap, over grid, between
 *            group and bands directly below and above it (\infty if group is
 *            whole spectrum)
 *
 * Returns TB_ERR_DEGENERATE (flux_out undefined) if that gap is below 1e-9 or a
 * link determinant is below 1e-10, i.e. when invariant is ill-defined or grid
 * is too coarse
 */
tb_status_t tb_berry_flux_grid(const tb_model_t *model, int first_band,
                               int n_bands, int n_k, double *flux_out,
                               double *min_gap);

/*
 * Chern number of band group (sum of tb_berry_flux_grid / 2 \pi), returned as a
 * double so distance to nearest integer shows whether grid is converged. 2D
 * models only. Same error conditions as tb_berry_flux_grid
 */
tb_status_t tb_chern_number(const tb_model_t *model, int first_band,
                            int n_bands, int n_k, double *chern_out);

/*
 * Wilson-loop (Berry / Zak) phase of band group along k1 at fixed k2 = k_other
 * (ignored for 1D models), in (-\pi, \pi]. For a 1D chain this is Zak phase.
 * Same degeneracy conditions as flux grid
 */
tb_status_t tb_wilson_phase(const tb_model_t *model, int first_band,
                            int n_bands, double k_other, int n_k,
                            double *phase_out);

/*
 * Winding of Wilson-loop phase as k2 runs once around zone:
 * (1/2 \pi) \sum \Delta\gamma, with `n_k2` steps (>= 4) in k2 and `n_k1` points
 * in each loop. For a gapped group it equals Chern number, giving a second,
 * independent way to compute it ("Wilson-loop flow" used to diagnose topology
 * from hybrid Wannier centers). Returns TB_ERR_DEGENERATE if phase jumps by
 * more than 0.9 \pi between consecutive k2 values (refine n_k2). 2D models only
 */
tb_status_t tb_wilson_winding(const tb_model_t *model, int first_band,
                              int n_bands, int n_k1, int n_k2,
                              double *winding_out);

/* Ready-made models (all return NULL on invalid / non-finite input)   */

// 1D chain, one orbital: E(k) = 2 t \cos k
tb_model_t *tb_model_chain(double t);

// Su-Schrieffer-Heeger chain, two orbitals: intracell hopping v, intercell w
// (H_01 = v + w \exp^{-i k}). Topological (Zak phase \pi) when |w| > |v|
tb_model_t *tb_model_ssh(double v, double w);

// 2D square lattice, one orbital: E = 2 t (\cos(k1) + \cos(k2))
tb_model_t *tb_model_square(double t);

// Honeycomb lattice (graphene), two orbitals (A, B), nearest-neighbor hopping
// t: H_AB = t (1 + \exp^{-i k1} + \exp^{-i k2}); Dirac points at K = (2 \pi/3,
// -2 \pi/3) and K' = -K
tb_model_t *tb_model_graphene(double t);

/*
 * Haldane model on honeycomb lattice: nearest-neighbor t1, complex
 * second-neighbor hopping t2 e^{\pm i phi} (+ on sublattice A, - on B, along
 * (n1, n2) = (1, 0), (-1, 1), (0, -1)), and a sublattice mass +M on A and -M on
 * B. A Chern insulator (lower band C = sign(t2 \sin(\phi))) when
 * |M| < 3 \sqrt{3} |t2 \sin(\phi)| and trivial otherwise. Haldane's analysis
 * assumes |t2| <= |t1|/3.
 */
tb_model_t *tb_model_haldane(double t1, double t2, double phi, double mass);

// Kagome lattice, three orbitals, nearest-neighbor hopping t. One exactly flat
// band at E = -2 t; other two touch it (and each other, at a Dirac point
// with E = t)
tb_model_t *tb_model_kagome(double t);

/*
 * Hofstadter model: square lattice in a uniform magnetic field with flux
 * p/q per plaquette (1 <= p < q), as a q-orbital magnetic unit cell. Built to
 * be identical to lattice_hofstadter_bloch() in lattice.h. Its q bands carry
 * Chern numbers given by TKNN Diophantine equation
 */
tb_model_t *tb_model_hofstadter(int p, int q, double t);

// Analytic references
double tb_chain_energy(double t, double k);
double tb_square_energy(double t, double k1, double k2);

// Upper graphene band |t| |1 + \exp^{-i k1} + \exp^{-i k2}| (lower band is
// its negative)
double tb_graphene_energy(double t, double k1, double k2);

/*
 * Gap of Haldane model at Dirac valley K (valley = +1,
 * k = (2 \pi/3, -2 \pi/3)) or K' (valley = -1):
 *   2 |M - valley 3 \sqrt{3} t2 \sin(\phi)|
 *
 * Returns NaN for any other valley value or non-finite input
 */
double tb_haldane_valley_gap(double t2, double phi, double mass, int valley);

/*
 * Chern number of lower Haldane band from phase diagram: sign(t2
 * \sin(\phi)) if |M| < 3 \sqrt{3} |t2 \sin(\phi)|, 0 if larger. Returns
 * TB_CHERN_UNDEFINED on a phase boundary (gapless), outside |t2| <= |t1|/3, for
 * t1 = 0, or for non-finite input
 */
int tb_haldane_chern(double t1, double t2, double phi, double mass);

/*
 * Zak phase of lower SSH band in this module's gauge: 0 if |v| > |w| (trivial),
 * \pi if |v| < |w| (topological). NaN at |v| == |w| (gapless) or for non-finite
 * input
 */
double tb_ssh_zak_phase(double v, double w);

// Energy of kagome flat band, -2 t
double tb_kagome_flat_band_energy(double t);

#endif // QMC_TIGHT_BINDING_H
