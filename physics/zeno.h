/*
 * Quantum Zeno effect for a driven two-level system
 *
 * A resonant drive with Rabi frequency Omega flips the qubit |0> -> |1> in
 * time \pi / Omega. Watching it suppresses the flip:
 *
 *  - Pulsed: N projective measurements of "still in |0>?" spaced T / N apart
 *    Each interval survives with probability \cos^2(\Omega T / (2N)), so total
 *    survival probability is cos^{2N}(Omega T / 2N) -> exp(-(Omega T)^2 / 4N)
 *    -> 1
 *  - Continuous: measurement dephases qubit at rate \Gamma
 *    Bloch equations dz/dt = -\Omega y, dy/dt = \Omega z - \Gamma y (z = P0 -
 *    P1) have slowest decay rate \Gamma / 2 for \Gamma <= 2 \Omega and (\Gamma
 *    - \sqrt(\Gamma^2 - 4 \Omega^2)) / 2 ~ \Omega^2 / \Gamma above it: decay is
 *    fastest at \Gamma = 2 \Omega and a stronger measurement slows it (Zeno)
 *    while a weaker one also does (no coupling to the apparatus)
 */
#ifndef ZENO_H
#define ZENO_H

/* \cos^{2N}(\omega T / 2N); NaN unless N >= 1, T >= 0 and omega finite */
double zeno_survival_pulsed(double omega, double t_total, int n_meas);

/*
 * Survival probability at any time 0 <= t <= T during pulsed experiment:
 * p^k \cos^2(\omega (t - k \tau) / 2) with \tau = T / N and k = floor(t / \tau)
 * measurements already passed. NaN on bad input
 */
double zeno_survival_pulsed_at(double omega, double t_total, int n_meas,
                               double t);

/*
 * Probability of still being in |0> at time t under continuous dephasing \gamma
 * (Bloch equations, RK4); NaN on bad input */
double zeno_survival_continuous(double omega, double gamma, double t);

/*
 * Slowest relaxation rate of Bloch equations (see above); NaN on bad input
 */
double zeno_continuous_rate(double omega, double gamma);

/*
 * Monte-Carlo trajectories of pulsed experiment: n_traj qubits, N measurements
 * alive[k] (k = 0 .. N) receives number of trajectories that passed all of
 * first k measurements. Returns 0, or -1 on bad input
 */
int zeno_sample_pulsed(double omega, double t_total, int n_meas, int n_traj,
                       unsigned seed, int *alive);

#endif
