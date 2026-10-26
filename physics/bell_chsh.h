/*
 * CHSH Bell test on a two-qubit Werner state
 *
 * \rho = V |\psi-><\psi-| + (1 - V) I / 4, singlet mixed with white noise
 * (visibility V). Alice and Bob measure spin along directions in x-z plane
 * at angles a and b, A(a) = \cos(a) \sigma_z + \sin(a) \sigma_x. Their
 * correlation is E(a, b) = Tr[\rho A(a) (x) A(b)] = -V \cos(a - b), computed
 * here from density matrix, not from closed form
 *
 * CHSH quantity S = E(a,b) - E(a,b') + E(a',b) + E(a',b'). Any local hidden
 * variable theory has |S| <= 2; quantum mechanics reaches 2 \sqrt{2} V
 * (Tsirelson), so Bell inequality is violated for V > 1 / \sqrt{2}
 * A concrete local model (a hidden angle, each side outputs sign of projection)
 * gives sawtooth correlation bell_local_correlation
 */

#ifndef BELL_CHSH_H
#define BELL_CHSH_H

/* E(a, b) from density matrix; NaN unless 0 <= V <= 1 and angles are finite */
double bell_correlation(double visibility, double a, double b);

/* S = E(a,b) - E(a,b2) + E(a2,b) + E(a2,b2); NaN on bad input */
double bell_chsh(double visibility, double a, double a2, double b, double b2);

/* Correlation of hidden-angle local model: -(1 - 2 d / \pi), d = |a - b| folded
 * into [0, \pi] */
double bell_local_correlation(double a, double b);

/* CHSH value of that local model for same settings */
double bell_local_chsh(double a, double a2, double b, double b2);

/* Largest |S| over 16 deterministic local strategies (always 2) */
double bell_deterministic_bound(void);

/* Largest |S| found by scanning a 4-angle grid with n_grid steps per angle over
 * [0, 2 \pi); NaN on bad input */
double bell_max_chsh(double visibility, int n_grid);

/* Simulate n_pairs measurements for each of four setting pairs (a,b), (a,b2),
 * (a2,b), (a2,b2) with a deterministic generator (seed)
 * outcome_sums[i] gets sum of products A*B (so E_est = outcome_sums[i] /
 * n_pairs)
 *
 * Returns estimate of S, or NaN on bad input */
double bell_sample(double visibility, double a, double a2, double b, double b2,
                   int n_pairs, unsigned seed, double outcome_sums[4]);

#endif
