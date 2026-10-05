#ifndef QMC_DIRAC_EVOLVE_H
#define QMC_DIRAC_EVOLVE_H

#include "../core/vector.h"

/*
 * Time evolution of 1D Dirac equation (two-component spinor), by split-operator
 * Fourier method
 *
 *   i \hbar d\psi/dt = [ c \sigma_x p + m c^2 \sigma_z + V(x) ] \psi,
 *   \psi = (upper, lower)^T
 *
 * Each step is a Strang splitting
 *   \psi(t + dt) = \exp^{-i V dt/(2\hbar)} . \exp^{-i H_0 dt/\hbar} . \exp^{-i
 * V dt/(2\hbar)} \psi(t) Free propagator is exact in momentum space, where
 * H_0(k) is a 2x2 matrix with H_0^2 = E_k^2 (E_k = \sqrt{(c \hbar k)^2 + (m
 * c^2)^2}): e^{-i H_0 dt/\hbar} = cos(E_k dt/\hbar) - i sin(E_k dt/\hbar)
 * H_0(k)/E_k so no Trotter error comes from kinetic/mass terms (they do not
 * commute with each other but are combined exactly), and only error is usual
 * O(dt^2) from V. scheme is unitary: total probability is conserved to
 * round-Off
 *
 * NOTE: Grid and boundaries. Arrays have N = 2^n points (N >= 4, radix-2 FFT)
 * with spacing dx. Point i sits at x_i = (i - N/2) dx. The FFT makes the grid
 * periodic, so a packet leaving one edge re-enters at the other: keep packets
 * away from the edges (or add absorption in the caller). A potential that is
 * not periodic (a step, say) has a discontinuity at the wrap-around point too;
 * this is harmless as long as nothing reaches it. Momentum grid is k_j = 2 \pi
 * j / (N dx) for j < N/2 and 2 \pi (j - N) / (N dx) otherwise, fft() uses
 * \exp^{-ikx}, so \exp^{i k0 x} has momentum +k0
 *
 * Spatial resolution matters for sharp potentials: waves inside a high step
 * have wave number up to \sqrt{(V0 - E)^2 - m^2 c^4}/(\hbar c), which must stay
 * well below \pi/dx
 */

/*
 * Advance (upper, lower) by `steps` steps of size dt under potential V
 * (length N, same units as energy; NULL means V = 0). Spinor is modified in
 * place
 *
 * Returns 0 on success, -1 on invalid input (NULL spinor, mismatched or
 * non-power-of-two length, dx <= 0, dt <= 0, steps < 1, m < 0, hbar <= 0, c <=
 * 0, non-finite V), -2 on allocation failure. On failure the spinor is
 * untouched
 */
int dirac_evolve_1d(cvector_t *upper, cvector_t *lower, const double *V,
                    double dx, double dt, int steps, double m, double hbar,
                    double c);

/*
 * Fill (upper, lower) with a Gaussian wavepacket centred at x0 with mean
 * wave number k0 and position spread \sigma (|\psi|^2 has standard deviation
 * \sigma), normalised so dx * sum(|upper|^2 + |lower|^2) = 1
 *
 * branch selects the spinor structure:
 *   +1: positive-energy states only. Gaussian * spinor of E = +E_k0 is
 *       projected exactly onto positive-energy subspace in momentum space, so
 *       packet moves with group velocity +c^2 \hbar k / E and does not
 *       tremble.
 *   -1: negative-energy states only (E = -E_k). Same construction. These
 *       "antiparticle-like" packets move with velocity -c^2 \hbar k / |E|,
 *       opposite to k0.
 *   0 : no projection, spinor (1, 0) (upper component only). Mixture of both
 *       energy signs with equal weight at k0 = 0, which produces Zitterbewegung
 *
 * Returns 0 on success, -1 on invalid input (as dirac_evolve_1d, plus \sigma <=
 * 0, branch outside {-1, 0, 1}, non-finite x0/k0), -2 on allocation failure
 */
int dirac_packet_1d(cvector_t *upper, cvector_t *lower, double dx, double x0,
                    double k0, double sigma, int branch, double m, double hbar,
                    double c);

/* Total probability dx * sum(|upper|^2 + |lower|^2); NaN on invalid input */
double dirac_norm_1d(const cvector_t *upper, const cvector_t *lower, double dx);

/* Mean position <x> = \int x (|upper|^2 + |lower|^2) dx / norm, with x_i =
 * (i - N/2) dx; NaN on invalid input or zero norm */
double dirac_position_1d(const cvector_t *upper, const cvector_t *lower,
                         double dx);

/*
 * Energy expectation value <H> = <\psi| c \sigma_x p + m c^2 \sigma_z + V
 * |\psi> / <\psi|\psi> (kinetic part evaluated in momentum space, so it is
 * exact for the sampled wavefunction). V may be NULL. It is conserved by
 * dirac_evolve_1d to O(dt^2) (exactly for V = 0). NaN on invalid input
 */
double dirac_energy_1d(const cvector_t *upper, const cvector_t *lower,
                       const double *V, double dx, double m, double hbar,
                       double c);

/*
 * Analytic transmission probability of a plane wave of energy E > m c^2,
 * incident from the left on a sharp potential step V(x) = V0 for x > 0 (0
 * for x < 0), from continuity of the spinor at x = 0:
 *   k = \sqrt{E^2 - m^2 c^4}/(\hbar c),  q = sgn(E') \sqrt{E'^2 - m^2
 * c^4}/(\hbar c),  E' = E - V0
 *   \kappa = (q/k) (E + m c^2) / (E' + m c^2),   T = 4 \kappa / (1 + \kappa)^2,
 * R = 1 - T (q carries the sign of E' so that the transmitted wave has positive
 * group velocity.) Three regimes: E' > m c^2          : ordinary transmission,
 * T -> 1 as V0 -> 0 |E'| < m c^2        : no propagating state in the step, T =
 * 0 E' < -m c^2         : Klein regime (V0 > E + m c^2). The transmitted wave
 * is a negative-energy state; T does NOT vanish as V0 grows (V0 -> infinity
 * gives \kappa -> \sqrt{(E + mc^2)/(E - mc^2)}), unlike non-relativistic
 * tunnelling which falls exponentially with barrier height Probability here is
 * the transmitted share of the incident probability current. A finite-width
 * packet reproduces it away from the thresholds V0 = E -+ m c^2, where T(E)
 * varies quickly
 *
 * Returns NaN if E <= m c^2 or any argument is non-finite / m < 0 / hbar <= 0
 * / c <= 0
 */
double dirac_step_transmission(double E, double V0, double m, double hbar,
                               double c);

#endif // QMC_DIRAC_EVOLVE_H
