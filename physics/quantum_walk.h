#ifndef QMC_QUANTUM_WALK_H
#define QMC_QUANTUM_WALK_H

#include "../core/vector.h"

/*
 * Discrete-time quantum walk on a line
 *
 * A walker with a two-state coin (up, down) at integer sites x = j - n/2
 * One step applies coin to every site, then shifts "up" one site to right and
 * "down" one site to left:
 *
 *   [a']   [ \cos(th)  \sin(th)] [a]       a_new(x) = a'(x - 1)
 *   [b'] = [ \sin(th) -\cos(th)] [b]       b_new(x) = b'(x + 1)
 *
 * th = \pi/4 is Hadamard walk. Map is unitary, so total probability is conserved
 * to round-off. Unlike classical random walk, whose width grows as sqrt(t),
 * quantum walk spreads ballistically: for a walker started with an unbiased
 * coin variance obeys
 *
 *     <x^2> -> (1 - \sin(th)) t^2          (t >> 1)
 *
 * and probability piles up near edges x = +/- t \cos(th) of light cone
 * ("fronts"), with a lower density in middle. Classical walk, with same coin
 * flip probability 1/2 per step, has <x^2> = t
 *
 * NOTE: Arrays have n sites and wrap around (periodic). Keep number of steps
 * below n/2 so probability from left edge has not reached right edge: walker
 * can move one site per step
 */

#define QWALK_COIN_UP 0   /* start in |up> */
#define QWALK_COIN_DOWN 1 /* start in |down> */
#define QWALK_COIN_SYM 2  /* start in (|up> + i |down>)/\sqrt(2), symmetric */

/* Put walker at centre site (n/2) with chosen coin state
 *
 * Returns 0 on success, -1 on NULL, mismatched or n < 4 lengths, or bad `coin`
 */
int qwalk_init(cvector_t *up, cvector_t *down, int coin);

/* Advance `steps` steps with coin angle \theta (radians).
 *
 * Returns 0 on success, -1 on invalid input (NULL, mismatched length,
 * non-finite theta, steps < 1), -2 on allocation failure. On failure state is
 * untouched */
int qwalk_step(cvector_t *up, cvector_t *down, double theta, int steps);

/* Total probability sum(|up|^2 + |down|^2); NaN on invalid input */
double qwalk_norm(const cvector_t *up, const cvector_t *down);

/* P(x_j) = |up_j|^2 + |down_j|^2 into prob[0..n-1]. Returns 0 or -1 */
int qwalk_probability(const cvector_t *up, const cvector_t *down, double *prob);

/* Mean position and variance (in sites, centre site = 0); NaN on invalid
 * input or zero norm */
double qwalk_mean(const cvector_t *up, const cvector_t *down);
double qwalk_variance(const cvector_t *up, const cvector_t *down);

/*
 * Classical unbiased random walk after t steps from centre: P(x) =
 * C(t, (t + x)/2) / 2^t for x of same parity as t. Fills prob[0..n-1]
 * (index j is site j - n/2)
 *
 * Returns 0, or -1 for NULL prob, n < 4, t < 0 or t >= n/2
 */
int qwalk_classical(double *prob, int n, int t);

/* Limit of <x^2>/t^2 for quantum walk, [1 - \sin(\theta)] (\theta in [0, pi/2]);
 * NaN for non-finite input */
double qwalk_asymptotic_variance(double theta);

#endif // QMC_QUANTUM_WALK_H
