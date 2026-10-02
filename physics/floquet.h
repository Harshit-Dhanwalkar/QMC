#ifndef QMC_FLOQUET_H
#define QMC_FLOQUET_H

#include "../core/matrix.h"
#include "../core/vector.h"

/*
 * Floquet theory for periodically driven quantum systems:
 *   i d\psi/dt = H(t) \psi,   H(t + T) = H(t),   T = 2 \pi / \omega.
 *
 * Floquet theorem: there are n solutions
 *   \psi_\alpha(t) = e^{-i \epsilon_\alpha t} u_\alpha(t),   u_\alpha(t + T) =
 *   u_\alpha(t),
 * with quasi-energies \epsilon_\alpha defined modulo \omega ("Brillouin zone"
 * of drive). Two independent solvers are provided:
 *
 *  1. Time-domain (floquet_solve_time): integrate one-period propagator
 *     U(T) for any callback H(t) and diagonalize it,
 *       U(T) u_\alpha(0) = e^{-i \epsilon_\alpha T} u_\alpha(0)
 *  2. Sambe space (floquet_solve_sambe): expand H(t) = \sum_m H_m e^{i m \omega
 *     t}, build time-independent extended ("Floquet") Hamiltonian on
 *     (harmonic) x (system) space
 *       [H_F]_{(m a),(n b)} = (H_{m-n})_{ab} + m \omega \delta_{mn}
 *       \delta_{ab}
 *     truncated to |m| <= K, and diagonalize it with Hermitian solver
 *
 * NOTE: Convention: quasi-energies are folded into [-\omega/2, \omega/2) and
 * returned in ascending order. Floquet modes are returned at t = 0 as  columns
 * of `modes`, ordered like quasi-energies. Natural units (\hbar = 1)
 *
 * NOTE: Time-domain solver: U(T) is propagated by classic RK4 on matrix
 * ODE dU/dt = -i H(t) U (error O(h^4), h = T / steps). Since U(T) is unitary
 * (hence normal), it is diagonalized through Hermitian matrix
 *   A = (U + U^\dagger)/2 + s (U - U^\dagger)/(2i)
 * its real and imaginary parts commute with each other, so they share
 * eigenvectors, and a generic s separates all eigenphases. Each eigenvalue is
 * recovered from Rayleigh quotient <v|U|v>. If residual max |U v - \lambda v|
 * exceeds 1e-8 for one s (an accidental near-degeneracy of A), several other
 * values of s are tried and result of lowest-residual attempt is kept; achieved
 * residual is always reported in `eigen_residual`
 *
 * NOTE: Sambe solver: extended spectrum is quasi-energy set shifted by all
 * integer multiples of \omega. Sorted ascending, exact spectrum is a periodic
 * sequence of n classes (one replica of each class per period \omega), so ANY n
 * consecutive eigenvalues contain each class exactly once. n eigenvalues in
 * middle of truncated spectrum are taken (indices K*n .. K*n + n - 1), where
 * truncation errors are smallest; this needs no heuristic about which replica
 * is "physical" and is safe at exact degeneracies and at resonances, where
 * eigenvectors straddle harmonic blocks
 * Mode u_\alpha(0) is sum of that eigenvector's harmonic blocks (a replica
 * shift m -> m + 1 only relabels blocks, so sum, and folded quasi-energy, are
 * replica independent)
 * `truncation_weight` is largest probability of any selected eigenvector in
 * two outermost harmonic blocks: it should be ~ 0; increase K until it is.
 *
 * Reference: Shirley, Phys. Rev. 138, B979 (1965); Sambe, Phys. Rev. A 7, 2203
 * (1973); Grifoni and Hanggi, Phys. Rep. 304, 229 (1998)
 */

// Fill the n x n matrix H_out (preallocated, overwritten) with H(t).
typedef void (*floquet_hamiltonian_fn)(double t, void *params,
                                       cmatrix_t *H_out);

typedef struct {
  int n;                    // system dimension
  double omega;             // drive angular frequency
  double *quasienergies;    // n values in [-\omega/2, \omega/2), ascending
  cmatrix_t *modes;         // n x n, column \alpha = u_\alpha(0)
  double unitarity_error;   // time-domain: max |U^\dagger U - I| (0 for Sambe)
  double eigen_residual;    // time-domain: max |U v - \lambda v| (0 for Sambe)
  double truncation_weight; // Sambe: largest weight in outermost blocks (0 for
                            // time-domain)
} floquet_result_t;

/*
 * Time-domain Floquet solver.
 *
 * n     : system dimension (>= 1)
 * H     : callback returning H(t); must be Hermitian for physical results
 * omega : drive angular frequency (> 0, finite)
 * steps : RK4 steps per period (>= 4); error scales as steps^-4
 *
 * Returns a heap-allocated result, or NULL on invalid input or allocation
 * failure. Free with floquet_result_free.
 */
floquet_result_t *floquet_solve_time(int n, floquet_hamiltonian_fn H,
                                     void *params, double omega, int steps);

/*
 * Sambe-space Floquet solver for H(t) = \sum_{m=-M}^{M} H_m e^{i m \omega t}.
 *
 * n         : system dimension (>= 1)
 * harmonics : array of 2M + 1 n x n matrices ordered m = -M, ..., 0, ..., +M
 *             (harmonics[M] is static part). Hermiticity of H(t) requires
 *             H_{-m} = H_m^\dagger
 * M         : highest harmonic present (>= 0)
 * omega     : drive angular frequency (> 0, finite)
 * K         : harmonic truncation, |m| <= K (K >= M + 1; larger is more
 *             accurate, extended matrix has n (2K + 1) rows)
 *
 * Returns a heap-allocated result, or NULL on invalid input or allocation
 * failure. Free with floquet_result_free
 */
floquet_result_t *floquet_solve_sambe(int n, cmatrix_t *const *harmonics, int M,
                                      double omega, int K);

void floquet_result_free(floquet_result_t *result);

/*
 * Stroboscopic evolution from a time-domain result: state after
 * `n_periods` full periods (t = n_periods * T),
 *   \psi(nT) = \sum_\alpha <u_\alpha|\psi_0> e^{-i \epsilon_\alpha n T}
 *              u_\alpha(0)
 * Exact for any integer n_periods (>= 0), at cost O(n^2) independent of
 * n_periods. psi0 must have length result->n. Returns a newly allocated vector,
 * or NULL on invalid input
 */
cvector_t *floquet_stroboscopic_state(const floquet_result_t *result,
                                      const cvector_t *psi0, long n_periods);

/*
 * Circular-drive two-level system (exactly solvable, no RWA):
 *   H(t) = (\omega_0/2) \sigma_z + (\Omega/2) (\cos(\omega t) \sigma_x +
 *                                              \sin(\omega t) \sigma_y).
 * In frame co-rotating with drive H is static, giving generalized
 * Rabi frequency \Omega_R = \sqrt{(\omega_0 - \omega)^2 + \Omega^2} and
 * quasi-energies \omega/2 \pm \Omega_R/2 (mod \omega). Returns folded
 * quasi-energy \omega/2 + sign * \Omega_R/2 into [-\omega/2, \omega/2); pass
 * sign = +1 or -1
 */
double floquet_circular_rabi_quasienergy(double omega0, double Omega,
                                         double omega, int sign);

/*
 * Cylindrical Bessel function J_0(x) (integral representation, trapezoid rule
 * with exponentially convergent accuracy ~1e-15 for |x| <= 50). Exposed because
 * coherent-destruction-of-tunneling prediction below needs it
 */
double floquet_bessel_j0(double x);

/*
 * High-frequency (\omega >> \Delta) quasi-energy splitting of periodically
 * driven double well / two-level system
 *   H(t) = (\Delta/2) \sigma_x + (A/2) \cos(\omega t) \sigma_z :
 *   splitting = \Delta |J_0(A/\omega)|   (Grossmann et al. 1991).
 * It vanishes at zeros of J_0 (coherent destruction of tunneling),
 * A/\omega = 2.4048, 5.5201, ... Corrections are O(\Delta^2/\omega)
 */
double floquet_cdt_splitting(double Delta, double A, double omega);

#endif // QMC_FLOQUET_H
