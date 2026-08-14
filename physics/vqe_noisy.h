#ifndef QMC_VQE_NOISY_H
#define QMC_VQE_NOISY_H

#include "../core/matrix.h"

/*
 *  NOTE: Noisy VQE: the same ansatz as vqe.c (vqe_prepare_ansatz: per layer,
 * RY(\theta) on every qubit then a linear CNOT chain), but simulated as an NxN
 * density matrix under single-qubit amplitude-damping (T1) and dephasing (T2)
 * noise applied after every gate, reusing lindblad.c's already-validated
 * channel builders and RK4 stepper.
 *
 * Noise model: after each single-qubit gate and after each CNOT in the ansatz,
 * the density matrix is evolved for a fixed duration `gate_time` under
 * lindblad_step_rk4 with H=0 (so only the dissipators act, not a unitary drift)
 * and L = one amplitude-damping + one dephasing operator per qubit (2 *
 * n_qubits operators total, at rates \gamma1/\gamma2). This is a
 * idle/gate-error noise model for NISQ-style simulation: everything remains
 * trace-preserving (Tr(\rho)=1 throughout) since each step is an exact GKSL
 * evolution, not an ad-hoc decay.
 *
 * \gamma1, \gamma2, gate_time are all in natural units (1/time and time,
 * \hbar=1); \gamma1=\gamma2=0 exactly reproduces the noiseless pure-state
 * result from vqe.c's vqe_energy.
 */

/*
 * Build the NxN (N = 2^n_qubits) noisy density matrix produced by running
 * ansatz with parameters theta under the noise model above.
 *
 * Returns NULL on invalid input (n_qubits<1, n_layers<1, !theta, \gamma1<0,
 * \gamma2<0, gate_time<0) or allocation failure.
 */
cmatrix_t *vqe_noisy_prepare_density(int n_qubits, int n_layers,
                                     const double *theta, double gamma1,
                                     double gamma2, double gate_time);

/*
 * Tr(H * \rho) for the noisy ansatz state at parameters \theta: noisy analogue
 * of vqe.c's vqe_energy.
 *
 * Returns 0.0 on invalid input (see vqe_noisy_prepare_density) or if H's
 * dimensions don't match 2^n_qubits.
 */
double vqe_noisy_energy(int n_qubits, int n_layers, const double *theta,
                        const cmatrix_t *H, double gamma1, double gamma2,
                        double gate_time);

/*
 *  NOTE: Zero-noise extrapolation (ZNE): evaluate vqe_noisy_energy at noise
 * scale factors c = 1, 2, ..., n_scales (\gamma1, \gamma2, and gate_time all
 * scaled by c; amplifying total decoherence proportionally, same qualitative
 * effect real-hardware ZNE achieves by stretching gate pulses or inserting
 * identity pairs), then fits a degree-(n_scales - 1) polynomial through those
 * n_scales points and extrapolates back to c=0.
 *
 * n_scales=2 gives the theoretical two-point linear extrapolation. n_scales=3
 * (default) fits a quadratic instead; empirically, for ansatz parameters
 * already at a VQE-optimized stationary point. A linear fit through such a
 * curve systematically overshoots past the true value rather than approaching
 * it.
 * n_scales must be >= 2.

 * Returns the extrapolated energy, and (if raw_c1 is non-NULL) writes
 * unmitigated c=1 energy to *raw_c1 for comparison.
 * Returns 0.0 (and leaves *raw_c1 untouched) on invalid input.
 */
double vqe_noisy_zne_energy(int n_qubits, int n_layers, const double *theta,
                            const cmatrix_t *H, double gamma1, double gamma2,
                            double gate_time, int n_scales, double *raw_c1);

#endif
