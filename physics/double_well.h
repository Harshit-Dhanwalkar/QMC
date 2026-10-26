/*
 * Tunnelling in a symmetric double well (ammonia-inversion physics)
 *
 *  H = -(1/2) d^2/dx^2 + V(x)
 * with \hbar = m = 1 on a
 *  Dirichlet grid x_i = -L + (i + 1) dx
 *  dx = 2L / (n + 1)
 *  V(x) = B (1 - x^2 / a^2)^2
 *
 * two wells at x = -+a separated by a barrier of height B
 *
 * Lowest eigenstates come in parity doublets (even, odd) whose splitting
 * \Delta = E1 - E0 shrinks exponentially with barrier
 * (WKB: \Delta ~ (\omega / \pi) \exp(-\int \sqrt(2 (V - E)) dx across barrier,
 * \omega = \sqrt(8 B / a^2) small-oscillation frequency of each well)
 * A particle prepared in one well, (\psi_0 - \psi_1) / \sqrt(2) (with sign
 * convention below), tunnels to other and back with P_left(t) = \cos^2(\Delta t
 * / 2): it oscillates with period 2 \pi / \Delta
 *
 * Solver is a tridiagonal Sturm bisection with parity-projected inverse
 * iteration, so near-degenerate doublets come out cleanly; potential must be
 * even
 */

#ifndef DOUBLE_WELL_H
#define DOUBLE_WELL_H

enum { DW_MAX_GRID = 4096, DW_MAX_STATES = 12 };

/* Fills v with V(x_i) = B (1 - x^2 / a^2)^2 on grid of n points over
 * [-half_width, half_width]; 0, or -1 on bad input */
int dw_potential(double barrier, double a, double half_width, int n, double *v);

/* Grid spacing for n points on [-L, L]. */
double dw_spacing(int n, double half_width);

/* Lowest n_states eigenvalues e[k] (ascending) and normalised eigenvectors
 * \psi[k * n + i] (\sum \psi^2 dx = 1) of an even potential v. Even states are
 * positive at their central-most maximum; odd states are positive on right
 * side. 0 on success, -1 on bad input (n <= DW_MAX_GRID,
 * n_states <= DW_MAX_STATES, v even to 1e-9 relative) */
int dw_eigen(const double *v, int n, double dx, int n_states, double *e,
             double *psi);

/* Overlaps c[k] = <\psi_k | g> of a normalised Gaussian g(x) ~ \exp(-(x -
 * center)^2 / (4 \sigma^2)) with eigenstates; *captured receives sum c^2 (1
 * when states span packet). 0 or -1 */
int dw_project(const double *psi, int n, double dx, int n_states,
               double half_width, double center, double sigma, double *c,
               double *captured);

/* \psi(t) = \sum_k c_k \exp(-i e_k t) \psi_k as (re, im) arrays of length n. 0
 * or -1 */
int dw_evolve(const double *psi, const double *e, const double *c, int n,
              int n_states, double t, double *re, double *im);

/* Probability to left of barrier (x < 0) of a state (re, im), normalised by
 * total probability */
double dw_prob_left(const double *re, const double *im, int n, double dx);

/* Mean position <x> for a state (re, im), normalised */
double dw_mean_x(const double *re, const double *im, int n, double dx,
                 double half_width);

#endif // DOUBLE_WELL_H
