#ifndef QMC_SSH_CHAIN_H
#define QMC_SSH_CHAIN_H

/*
 * Su-Schrieffer-Heeger (SSH) chain: the simplest topological insulator
 *
 * A finite chain of n cells, each with an A and a B site (2n sites in all, site
 * 2c = A of cell c, site 2c + 1 = B of cell c). Intra-cell hopping v joins
 * A_c and B_c, inter-cell hopping w joins B_c and A_(c+1):
 *
 *   H = sum_i eps_i |i><i| + sum_i bond_i (|i><i+1| + h.c.)
 *
 * with bond_(2c) = v and bond_(2c+1) = w for the clean chain. The bulk bands
 * are E(k) = +/- |v + w e^{ik}| = +/- sqrt(v^2 + w^2 + 2 v w cos k), with a gap
 * 2|v - w| that closes at |v| = |w|. The winding number of h(k) = v + w e^{ik}
 * around the origin is 1 for |w| > |v| (topological) and 0 for |w| < |v|
 * (trivial). In the topological phase a finite chain has two edge states
 * with energies ~ +/- (v/w)^n (exponentially close to zero, inside the gap):
 * one on the A sublattice at the left end with amplitude ~ (-v/w)^c on cell c,
 * its partner on the B sublattice at the right end. Their decay length is
 * xi = 1 / ln(w/v) cells.
 *
 * Chiral (sublattice) symmetry, eps = 0 and hoppings only between A and B,
 * makes the spectrum symmetric under E -> -E and pins the edge states at zero
 * energy up to the exponentially small left-right splitting: random hopping
 * does not move them. A staggered potential eps = (+m, -m, ...) breaks the
 * symmetry and shifts them to +/- m.
 */

/* Bond strengths bond[0..2n-2] of the chain: v on intra-cell bonds and w on
 * inter-cell bonds, each multiplied by (1 + noise u) with u uniform in
 * [-1, 1] from a deterministic generator seeded by `seed` (noise = 0 gives the
 * clean chain). Random hopping keeps the chain chiral symmetric. Returns 0, or
 * -1 for NULL bond, n_cells < 1, non-finite v or w, or noise outside
 * [0, 1). */
int ssh_chain_bonds(double *bond, int n_cells, double v, double w, double noise,
                    unsigned long long seed);

/* On-site energies eps[0..2n-1] = +m on A sites and -m on B sites. Returns 0,
 * or -1 for NULL eps, n_cells < 1 or non-finite m. */
int ssh_chain_onsite(double *eps, int n_cells, double m);

/*
 * Diagonalise the chain. bond has 2n-1 entries; eps has 2n entries or is NULL
 * (all zero). Writes the 2n eigenvalues in ascending order to energies and the
 * eigenvectors to vectors[s * 2n + i] = amplitude of state s on site i (real,
 * normalised); vectors may be NULL. Returns 0, or -1 on invalid input (NULL
 * bond or energies, n_cells < 1, non-finite entries), -2 on allocation
 * failure.
 */
int ssh_chain_solve(int n_cells, const double *bond, const double *eps,
                    double *energies, double *vectors);

/* Upper bulk band sqrt(v^2 + w^2 + 2 v w cos k); NaN for non-finite input. */
double ssh_bulk_energy(double v, double w, double k);

/* Winding number of v + w e^{ik} around the origin (1 if |w| > |v|, 0 if
 * |w| < |v|), computed by unwrapping the phase over a fine k grid. Returns -1
 * at the gap closing |v| = |w| and for non-finite input. */
int ssh_winding_number(double v, double w);

/* Zak phase of the lower band, |gamma| in [0, pi], from a discretised Berry
 * phase (0 for the trivial and pi for the topological phase). NaN for
 * non-finite input or v = w = 0. */
double ssh_zak_phase(double v, double w);

/* Edge-state decay length 1 / ln(|w|/|v|) in cells for |w| > |v| > 0; NaN
 * otherwise. */
double ssh_edge_decay_length(double v, double w);

#endif // QMC_SSH_CHAIN_H
