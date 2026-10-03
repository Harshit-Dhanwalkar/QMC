#ifndef QMC_TD_PERTURBATION_H
#define QMC_TD_PERTURBATION_H

#include "../core/complex.h"
#include "../core/matrix.h"
#include "../core/vector.h"

/*
 * Time-dependent perturbation theory (Dyson series)
 *
 * H(t) = H_0 + V(t), with H_0 diagonal in working basis
 * H_0 |n> = E_n |n>. Natural units (\hbar = 1)
 *
 * NOTE: Interaction picture. With Schrodinger amplitudes
 * a_n(t) = <n|\psi(t)> and interaction-picture amplitudes
 *   - c_n(t) = e^{i E_n t} a_n(t),
 * Schrodinger equation becomes
 *   - dc_f/dt = -i \sum_m e^{i (E_f - E_m) t} V_{fm}(t) c_m(t),
 * and iterating from c(t_0) = \psi_0 gives Dyson series
 *   - c = c^{(0)} + c^{(1)} + c^{(2)} + ...,   c^{(0)} = \psi_0
 *   - dc^{(k)}_f/dt = -i \sum_m e^{i (E_f - E_m) t} V_{fm}(t) c^{(k-1)}_m(t)
 *
 * Where
 *   - c^{(k)} is exactly k-th order in V (it scales as \epsilon^k when V ->
 *     \epsilon V)
 *   First two orders are theoritical results
 *   - c_f^{(1)}(t) = -i \int_{t_0}^t dt' \exp^{i \omega_{fi} t'} V_{fi}(t')
 *   - c_f^{(2)}(t) = -\sum_m \int_{t_0}^t dt_1 \int_{t_0}^{t_1} dt_2
 *                     \exp^{i \omega_{fm} t_1} V_{fm}(t_1)
 *                     \exp^{i \omega_{mi} t_2} V_{mi}(t_2)
 *
 * NOTE: All orders 1..K are integrated together as one coupled hierarchy of
 * linear ODEs with classic RK4 (error O(h^4), h = (t1 - t0)/steps), rather than
 * by nested quadratures, so cost is O(steps * K * n^2) for any K. Each RK4 step
 * calls perturbation callback at t, t + h/2 and t + h
 *
 * NOTE: Truncating series at order K leaves an error O(\epsilon^{K+1}):
 * truncated state is not exactly normalized, and tdpt_norm reports by how much.
 * Series converges only while perturbation is "small" over time window (|V| t
 * << 1 for resonant coupling, or any |V|/|\omega| << 1 off resonance); near
 * resonance at long times use a non-perturbative method
 *
 * Probabilities are picture independent: |c_n|^2 = |a_n|^2
 */

// Fill n x n matrix V_out (preallocated and zeroed before every call) with
// V(t) in eigenbasis of H_0. Must be Hermitian for physical results
typedef void (*tdpt_perturbation_fn)(double t, void *params, cmatrix_t *V_out);

typedef struct {
  int n;             // number of H_0 eigenstates
  int max_order;     // highest order computed (>= 0)
  double t0, t1;     // integration window; amplitudes are at t1
  double *energies;  // copy of n H_0 energies
  cvector_t *orders; // max_order + 1 vectors held by value; orders[k] =
                     // c^{(k)}(t1). Owned by result: never cvector_free an
                     // element, only call tdpt_result_free
} tdpt_result_t;

#define TDPT_MAX_DIM 1024
#define TDPT_MAX_ORDER 12

/*
 * Compute Dyson-series amplitudes c^{(0..max_order)}(t1)
 *
 * n         : number of states (1 .. TDPT_MAX_DIM)
 * energies  : H_0 eigenvalues E_n, length n (finite)
 * V         : perturbation callback, V(t) in H_0 eigenbasis
 * psi0      : initial state at t0 (length n); need not be a basis state
 * t0, t1    : integration window, t1 > t0, both finite
 * steps     : RK4 steps over window (>= 4)
 * max_order : highest order (0 .. TDPT_MAX_ORDER); 0 just returns psi0
 *
 * Returns a heap-allocated result, or NULL on invalid input, non-finite output,
 * or allocation failure
 */
tdpt_result_t *tdpt_dyson(int n, const double *energies, tdpt_perturbation_fn V,
                          void *params, const cvector_t *psi0, double t0,
                          double t1, int steps, int max_order);

void tdpt_result_free(tdpt_result_t *result);

/*
 * Interaction-picture amplitude summed through `order`:
 *   \sum_{k=0}^{order} c^{(k)}(t1), for 0 <= order <= max_order
 *
 * Returns a new vector, or NULL on invalid input
 */
cvector_t *tdpt_total_amplitude(const tdpt_result_t *result, int order);

/*
 * Same sum converted to Schrodinger picture:
 *   a_n(t1) = \exp^{-i E_n t1} \sum_{k <= order} c_n^{(k)}(t1)
 *
 * Returns a new vector, or NULL on invalid input
 */
cvector_t *tdpt_schrodinger_amplitude(const tdpt_result_t *result, int order);

/*
 * Occupation probability |\sum_{k <= order} c_state^{(k)}|^2 of basis state
 * `state`
 *
 * Returns NaN on invalid input (0 is a legitimate probability)
 */
double tdpt_probability(const tdpt_result_t *result, int order, int state);

/*
 * Norm \sum_n |\sum_{k <= order} c_n^{(k)}|^2 of truncated series
 * Equals initial norm only up to O(\epsilon^{order + 1}); deviation is a direct
 * measure of truncation error
 *
 * Returns NaN on invalid input
 */
double tdpt_norm(const tdpt_result_t *result, int order);

/*
 * Closed-form building blocks for first-order theory (ground truth for tests
 * and for quick estimates)
 */

/*
 * \int_0^t \exp^{i \Delta s} ds = (\exp^{i \Delta t} - 1) / (i \Delta), with
 * \Delta -> 0 limit (= t) handled by a series so result is continuous
 */
complex_t tdpt_phase_integral(double delta, double t);

/*
 * Constant perturbation V_fi switched on at t = 0:
 *   c_f^{(1)}(t) = -i V_fi \int_0^t \exp^{i \omega_{fi} s} ds
 * with \omega_{fi} = E_f - E_i, giving:
 *   P = 4 |V_fi|^2 \sin^2(\omega t / 2) / \omega^2
 */
complex_t tdpt_constant_first_order(complex_t Vfi, double omega_fi, double t);

/*
 * Harmonic perturbation V(t) = V cos(\Omega t) switched on at t = 0 (no RWA):
 *   c_f^{(1)}(t) = -i (V_fi / 2) [
 *                               \int_0^t \exp^{i (\omega_{fi} + \Omega) s} ds
 *                             + \int_0^t \exp^{i (\omega_{fi} - \Omega) s} ds
 *                                ]
 */
complex_t tdpt_harmonic_first_order(complex_t Vfi, double omega_fi,
                                    double Omega, double t);

/*
 * Gaussian pulse V(t) = V_0 \exp(-t^2 / (2 \tau^2)) acting over all time
 * (t: -inf -> +inf), final amplitude
 *   c_f^{(1)}(\infty) = -i V_fi \tau \sqrt{2 \pi}
 *                        * \exp^{-\omega_{fi}^2 \tau^2 / 2
 * adiabatic-suppression law: transitions with \omega \tau >> 1 are
 * exponentially unlikely
 *
 * Returns NaN for tau <= 0 or non-finite input
 */
complex_t tdpt_gaussian_first_order(complex_t Vfi, double omega_fi, double tau);

#endif // QMC_TD_PERTURBATION_H
