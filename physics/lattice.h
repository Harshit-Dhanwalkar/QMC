#ifndef QMC_LATTICE_H
#define QMC_LATTICE_H

#include "../core/matrix.h"
#include "../core/vector.h"
#include <stdint.h>

/*
 * Tight-binding lattice models (condensed matter): nearest-neighbor hopping
 * Hamiltonians on 1D/2D lattices, in the site (real-space) basis.
 *  H = \epsilon0 * \sum_i n_i - t * \sum_{<i,j>} (c_i^\dagger c_j + h.c.)
 *
 * With <i,j> summing over nearest-neighbor bonds.
 * NOTE: All Hamiltonians here are dense, real-valued (stored as complex_t with
 * zero imaginary part, for uniform use with this project's existing complex
 * eigensolvers), Hermitian matrices in site basis. Bloch's theorem / Fourier
 * analysis gives exact analytic eigenvalues for clean (disorder-free) cases,
 * which the lattice_*_analytic functions return for direct comparison against
 * numerical diagonalization.
 */

typedef enum {
  LATTICE_OPEN = 0,    /* open chain/boundary: no wraparound bond */
  LATTICE_PERIODIC = 1 /* periodic (ring/torus): wraparound bond included */
} lattice_bc_t;

/*
 * 1D tight-binding chain, n_sites sites, on-site energy \epsilon0, hopping
 * amplitude t.
 *
 * Returns an n_sites x n_sites dense Hermitian cmatrix_t, or NULL if n_sites
 * < 1.
 */
cmatrix_t *lattice_build_1d_chain(int n_sites, double epsilon0, double t,
                                  lattice_bc_t bc);

/*
 * Exact analytic eigenvalues of lattice_build_1d_chain (Bloch's theorem for
 * LATTICE_PERIODIC, standing-wave quantization for LATTICE_OPEN),
 * In ascending order into E_out:
 *   periodic: E_n = \epsilon0 - 2 * t * \cos(2 * \pi * n / n_sites),
 *                     n=0..n_sites-1
 *   open:     E_n = \epsilon0 - 2 * t * \cos(n * \pi / (n_sites + 1)),
 *                    n=1..n_sites
 */
void lattice_1d_chain_analytic(int n_sites, double epsilon0, double t,
                               lattice_bc_t bc, double *E_out);

/*
 * 2D square-lattice tight-binding, nx x ny sites (site (ix,iy) at row-major
 * index ix * ny + iy), same on-site epsilon0 and hopping t in both directions,
 * same boundary condition applied independently along x and y.
 *
 * Returns an (nx * ny) x (nx * ny) dense Hermitian cmatrix_t, or NULL if nx<1
 * or ny<1.
 */
cmatrix_t *lattice_build_2d_square(int nx, int ny, double epsilon0, double t,
                                   lattice_bc_t bc);

/*
 * Exact analytic eigenvalues of lattice_build_2d_square (2D lattice is exactly
 * separable into independent x- and y-direction 1D problems, since there are no
 * diagonal/cross hopping terms), ascending order, into E_out:
 *   E_{m,n} = \epsilon0 - 2 * t * [f(m,nx) + f(n,ny)]
 *  Where
 *   f = the periodic or open 1D dispersion argument from
 * lattice_1d_chain_analytic, minus its own \epsilon0 term
 */
void lattice_2d_square_analytic(int nx, int ny, double epsilon0, double t,
                                lattice_bc_t bc, double *E_out);

/*
 * 1D Anderson model: open-boundary tight-binding chain with uncorrelated
 * on-site disorder, on-site energy at site i drawn uniformly from
 * [-disorder_W/2, +disorder_W/2] (fixed by `seed` for reproducibility), hopping
 * amplitude t (uniform, no disorder in bonds). The paradigm model for Anderson
 * localization: in 1D, arbitrarily weak disorder (disorder_W > 0) is known to
 * localize every eigenstate, with a disorder-dependent localization length.
 *
 * Returns an n_sites x n_sites dense Hermitian cmatrix_t, or NULL if n_sites <
 * 1 or disorder_W < 0.
 */
cmatrix_t *lattice_build_anderson_1d(int n_sites, double t, double disorder_W,
                                     uint64_t seed);

/*
 * Inverse participation ratio IPR = \sum_i |\psi_i|^4 / (\sum_i |\psi_i|^2)^2,
 * standard localization diagnostic: IPR ~ 1/n_sites for a fully delocalized
 * (extended) state spread evenly over every site, IPR -> 1 for state localized
 * on single site. Qualitative signature used to distinguish extended from
 * localized states. Normalizes internally.
 *
 * Returns 0.0 if \psi is NULL or has zero norm.
 */
double lattice_ipr(const cvector_t *psi);

/*
 * Su-Schrieffer-Heeger (SSH) model: a dimerized 1D chain of n_cells unit cells
 * (2 * n_cells sites total), alternating intra-cell hopping t1 and inter-cell
 * hopping t2 (site 2c, 2c+1 are cell c's two sublattice sites; bond (2c,2c+1)
 * has amplitude -t1, bond (2c+1, 2c+2) has amplitude -t2). On-site energy is 0
 * for both sublattices (SSH convention; chiral/sublattice symmetry is makes
 * topological phase's near-zero-energy edge states exact zero modes in
 * infinite-chain limit).
 *
 * NOTE: Topological phase (t2 > t1, open boundary): two near-zero-energy states
 * appear, exponentially localized at two chain ends, signature of topological
 * edge state.
 *
 * Trivial phase (t1 > t2): no near-zero states, no edge localization.
 *
 * Where
 *  bc = LATTICE_PERIODIC closes chain into a ring (bond (2*n_cells-1, 0) with
 *       amplitude -t2), removing physical edges and therefore the edge states
 *       regardless of t1 vs t2 (topological invariant then only shows up in
 *       bulk band structure/Zak phase, not as edge states).
 *
 * Returns a (2 * n_cells) x (2 * n_cells) dense Hermitian cmatrix_t, or NULL if
 * n_cells < 1.
 */
cmatrix_t *lattice_build_ssh(int n_cells, double t1, double t2,
                             lattice_bc_t bc);

/*
 * Landau levels on a 2D square lattice via Peierls substitution (Hofstadter /
 * Harper model): a uniform perpendicular magnetic field is introduced by
 * multiplying hopping amplitudes by phase factors:
 *   \exp(i * e/\hbar * \integral of A.dl along the bond),
 * rather than adding any new on-site term.
 *
 * Landau gauge A = (0, B * x, 0): bonds along x are unaffected (amplitude -t,
 * real), bonds along y acquire a position-dependent phase set by site's x
 * index:
 *   H_{(ix, iy),(ix, iy + 1)} = -t * \exp(i * 2 * \pi * \alpha * ix)
 *   h.c. for the reverse bond
 *
 * Where
 *  \alpha = flux through one plaquette / flux quantum = B * a^2 / (2 * \pi) in
 * natural units (\hbar = e = 1, lattice constant a).
 * NOTE: \alpha is the Hofstadter-model flux parameter (\alpha = p/q rational
 * reproduces the Hofstadter butterfly at that flux fraction; \alpha to be
 * rational since it builds a finite open lattice, not a magnetic-unit-cell
 * Bloch Hamiltonian).
 *
 * NOTE: On boundary conditions: Landau gauge above breaks translational
 * symmetry along x (y-hopping phase depends explicitly on ix), so a naive
 * periodic wraparound bond in x would not be gauge-consistent on a finite
 * lattice; only bc_y controls periodicity along y (wraparound y-bonds carry
 * same phase rule, which is consistent). Therefore x is always open boundary.
 *
 * Returns an (nx * ny) x (nx * ny) dense Hermitian cmatrix_t, or NULL if nx<1
 * or ny<1.
 */
cmatrix_t *lattice_build_2d_square_magnetic(int nx, int ny, double epsilon0,
                                            double t, double alpha,
                                            lattice_bc_t bc_y);

/*
 * Continuum-limit (weak-field, alpha << 1) approximation to n-th Landau level
 * of lattice_build_2d_square_magnetic, measured from the zero-field band bottom
 * \epsilon_0 - 4 * t:
 *  E_n ~= \epsilon_0 - 4 * t + \omega_c * (n + 1/2),
 *
 * Where:
 *   \omega_c = 4 * \pi * t * \alpha
 *
 * NOTE: Derivation: expanding tight-binding dispersion near k=0 gives an
 * isotropic effective mass m* = 1 / (2 * t * a^2), the Peierls phase on y-bonds
 * expanded near a bond's x=0 gives a harmonic confining potential in x with
 * frequency :
 *  \omega_c = 2 * t * a^2 * B = 4 * pi * t * \alpha
 * (using \alpha = B * a^2 / (2 * \pi))
 * This is 2D magneto-tight-binding effective-mass result (References :
 * Hofstadter 1976; MacDonald 1983); each level is macroscopically degenerate in
 * \alpha<<1 limit (one state per flux quantum through the sample).
 *
 * Deviations at larger alpha come from lattice curvature (dispersion is not
 * perfectly parabolic) and, on a finite lattice, from edge effects once the
 * magnetic length ~1 / \sqrt(\alpha) becomes comparable to nx or ny.
 */
double lattice_landau_level_energy(int n, double epsilon0, double t,
                                   double alpha);

/*
 * Hofstadter model: 2D square-lattice tight-binding electron in a uniform
 * magnetic field with rational flux \alpha = p/q (p, q coprime integers; flux
 * quantum per plaquette). Unlike lattice_build_2d_square_magnetic (finite
 * real-space lattice, open/periodic boundary in y, used for weak-field
 * Landau-level limit above), this works directly in thermodynamic-limit
 * momentum-space (Bloch) picture: magnetic unit cell is q sites wide (Landau
 * gauge, phase carried on x-bonds as a function of site's row within cell,
 * period q in y), giving a q x q Bloch Hamiltonian H(kx,ky) whose q eigenvalues
 * are the q "Hofstadter butterfly" subbands at that (kx,ky). This is
 * momentum-space object the Chern-number calculation below integrates Berry
 * curvature over.
 *
 * Returns a newly-allocated q x q complex Hermitian matrix (NULL on invalid
 * input: p<1, q<2, or p>=q i.e. not a flux in (0,1)); caller owns it
 * (cmatrix_free).
 */
cmatrix_t *lattice_hofstadter_bloch(double kx, double ky, int p, int q,
                                    double t);

/*
 * Chern number of each of the q Hofstadter subbands, (Reference:
 * Fukui-Hatsugai-Suzuki (2005)) via discretized/gauge-invariant lattice
 * Berry-curvature method: sample H(kx,ky) on an n_k x n_k grid over (2 * \pi x
 * 2 * \pi) magnetic Brillouin zone, form the U(1) link variable U_{\mu}(k) =
 * <n,k|n,k + \mu> / |<n,k|n,k + \mu>| between neighboring grid points'
 * eigenvectors for each band n, and \sum arg(plaquette product of 4 links) over
 * every plaquette : this is gauge-invariant (no smooth-gauge choice needed
 * across the BZ, unlike a naive derivative-of-Berry-connection approach) and
 * gives an exactly-integer-quantized result in the n_k->\infty limit for any
 * gapped band.
 *
 * Writes q values into `chern_out` (caller-allocated, length q): each
 * should be within numerical-grid-resolution of an integer for any p/q
 * where subband n is gapped from its neighbors.
 *
 * NOTE: at band-touching points (e.g. p/q=1/2, the "pi-flux" case, where the 2
 * subbands touch at Dirac points) individual-band Chern numbers are not
 * well-defined by the usual TKNN argument - this function will still return a
 * number (bands near-touching but not exactly, on a finite grid, given some
 * residual curvature), but it need not be well-quantized in that case; the sum
 * over all q bands is always confined to a strong sum-rule (=0 exactly, a
 * useful self-consistency check) regardless of gap closures, since it is Chern
 * number of full Bloch bundle across compensating regions of curvature.
 *
 * Returns 1 on success, 0 on invalid input (p<1, q<2, p>=q, n_k<2, or
 * chern_out NULL).
 * NOTE: For even q, the Hofstadter spectrum can have exact band touchings
 * between adjacent subbands at isolated (kx,ky) points (Dirac points),
 * independent of numerical resolution : e.g. p/q=1/2 (all q=2 bands touch, both
 * individual Chern numbers 0) and p/q=1/4 (the middle two of four bands touch,
 * min gap exactly 0 across the whole BZ). At such a touching, the FHS
 * link-variable construction below is only meaningful away from the touching
 * point, and this function's per-band output can become numerically unstable /
 * grid-resolution-dependent right at touching bands specifically (verified: q=3
 * and q=5, both fully gapped for every kx,ky, give grid-independent integer
 * results from n_k=16 through n_k=100+; q=4's touching pair does not). The sum
 * over all q bands together remains exactly 0 regardless (full q-band bundle is
 * trivial). */
int lattice_hofstadter_chern_numbers(int p, int q, double t, int n_k,
                                     double *chern_out);

#endif
