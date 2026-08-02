#ifndef QMC_VQE_H
#define QMC_VQE_H

#include "../core/matrix.h"
#include "../core/vector.h"
#include <stdint.h>

/*
 * Variational Quantum Eigensolver (VQE): hybrid quantum-classical ground state
 * solver
 *
 * NOTE: Ansatz: hardware-efficient ansatz, n_layers of [RY(\theta) on every
 * qubit, then a linear chain of CNOTs for entanglement (qubit i controls qubit
 * i+1, i=0..n_qubits-2)]. Total parameters = n_qubits * n_layers.
 * Refrence : VQE ansatz (see e.g. Kandala et al. 2017).
 *
 * NOTE: Optimizer: coordinate descent, repeatedly sweep through every
 * parameter, minimizing each one via golden_section_minimize over a window
 * [\theta_i - window, \theta_i + window] around its current value with others
 * held fixed, for n_sweeps sweeps.
 *
 * NOTE: All Hamiltonians here are represented as dense Hermitian cmatrix_t of
 * size 2^n_qubits. There is no general Pauli-string/tensor-product operator
 * builder yet; only a transverse-field Ising model builder (vqe_build_tfim) is
 * present at this current stage. Fine for small qubit counts.
 * TODO: Pauli-string layer needed before scaling to anything larger or maps
 * molecular (Jordan-Wigner) Hamiltonians.
 */

typedef struct {
  double energy;     // <\psi(\theta_opt)|H|\psi(\theta_opt)>, final
  double *theta_opt; // n_qubits*n_layers optimized parameters
  int n_params;
} vqe_result_t;

/*
 * Prepare ansatz state |\psi(\theta)> = U(\theta)|00...0> for given parameters
 * (length n_qubits*n_layers, layer-major: \theta[layer * n_qubits + qubit]).
 *
 * Returns a allocated cvector_t, or NULL on invalid input.
 */
cvector_t *vqe_prepare_ansatz(int n_qubits, int n_layers, const double *theta);

/*
 * Exact expectation value Re(<\psi|H|\psi>) for a normalized state \psi and
 * Hermitian H (where H must be 2^n_qubits square, matching \psi's dimension).
 * Returns 0.0 on invalid input.
 */
double vqe_expectation(const cvector_t *psi, const cmatrix_t *H);

/*
 * Prepare ansatz for `theta` and return its energy expectation against H in one
 * call (frees intermediate state internally).
 */
double vqe_energy(int n_qubits, int n_layers, const double *theta,
                  const cmatrix_t *H);

/*
 * Run VQE: random-initialize \theta (fixed by `seed`, uniform in [-\pi, \pi]),
 * then coordinate-descend for n_sweeps sweeps. H must be a 2^n_qubits x
 * 2^n_qubits Hermitian cmatrix_t.
 *
 * Returns a vqe_result_t with final energy and optimized parameters. On invalid
 * input, returns a zeroed result with `theta_opt == NULL`.
 */
vqe_result_t vqe_run(int n_qubits, int n_layers, const cmatrix_t *H,
                     int n_sweeps, double window, uint64_t seed);

/*
 * Build dense (2^n_qubits x 2^n_qubits) Hamiltonian of 1D transverse-field
 * Ising model on an open chain of n_qubits sites:
 *   H = -J * \sum_{i = 0}^{n - 2} Z_i Z_{i + 1}  -  h * sum_{i = 0}^{n - 1} X_i
 *
 * Returns NULL if n_qubits < 1.
 */
cmatrix_t *vqe_build_tfim(int n_qubits, double J, double h);

#endif
