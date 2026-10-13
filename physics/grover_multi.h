/*
 * Grover search with several marked items, on full state vector
 *
 * N = 2^n items, M of them marked. One Grover iteration is oracle (flip sign of
 * every marked amplitude) followed by diffusion operator 2|s><s| - 1 (inversion
 * about mean). State stays in plane spanned by uniform marked and unmarked
 * superpositions and rotates by 2 \theta per iteration, \sin(\theta) = \sqrt(M
 * / N), so success probability after k iterations is \sin^2((2k + 1) \theta)
 * and best k is near \pi / (4 \theta) - 1/2
 *
 * State is real: amplitudes are stored as doubles
 */
#ifndef GROVER_MULTI_H
#define GROVER_MULTI_H

#define GROVER_MAX_QUBITS 12

/* amp[0 .. 2^n - 1] = uniform superposition 1 / \sqrt(N). 0 on success, -1 on
 * bad input */
int grover_multi_init(double *amp, int n_qubits);

/* One iteration; marked[i] != 0 flags a marked item. 0 on success, -1 on bad
 * input */
int grover_multi_step(double *amp, int n_qubits, const unsigned char *marked);

/* Total probability on marked items; NaN on bad input */
double grover_multi_success(const double *amp, int n_qubits,
                            const unsigned char *marked);

/* Rotation angle \theta = \asin(\sqrt(M / N)); NaN unless 1 <= M <= N (M = N
 * gives \pi / 2) */
double grover_multi_theta(int n_qubits, int n_marked);

/* Closed-form success probability \sin^2((2k + 1) \theta) */
double grover_multi_success_exact(int n_qubits, int n_marked, int k);

/* Iteration count with highest success probability, round(\pi / (4 \theta) -
 * 1/2); -1 on bad input */
int grover_multi_optimal_iterations(int n_qubits, int n_marked);

#endif
