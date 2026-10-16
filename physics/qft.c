/*
 * Quantum Fourier transform on a subset of qubits of a larger register
 *
 * Acts on a chosen set of n_q qubits inside an n_qubits-wide state vector,
 * with no requirement that qubits be contiguous or in any particular
 * order. Standard n_q-qubit QFT
 *
 *     |j> -> (1 / \sqrt{2^{n_q}}) \sum_{k=0}^{2^{n_q} - 1}
 *              \exp(2 \pi i j k / 2^{n_q}) |k>
 *
 * is decomposed as usual sandwich of Hadamards and controlled-phase gates,
 * followed by a bit-reversal permutation. Controlled-phase gate R_m between
 * qubit j (target) and qubit k (control), k > j in chosen ordering, has phase
 * 2\pi / 2^{k-j+1}; this makes gate count O(n_q^2), standard decomposition
 *
 * Inverse QFT is implemented as exact inverse gate sequence: reverse swaps,
 * then for each j in decreasing order apply controlled-phase gates with
 * opposite sign followed by an H. This avoids assuming anything about numerical
 * accuracy of forward transform and is what qpe.c relies on for phase-kickback
 * readout
 *
 * "Qubit j is most-significant" convention used by qpe.c is that first entry of
 * `qubits[]` is bit j = n_q - 1 of Fourier basis integer. This is what makes
 * QPE output a big-endian binary fraction of phase with no extra reversal step
 */
#include "qft.h"
#include "../core/vector.h"
#include "qubits.h"
#include <math.h>

void qft_apply(cvector_t *psi, int n_qubits, const int *qubits, int n_q) {
  if (!psi || !qubits || n_q < 1) {
    return;
  }

  for (int j = 0; j < n_q; j++) {
    qstate_apply_gate1(psi, n_qubits, qubits[j], hadamard_gate);

    for (int k = j + 1; k < n_q; k++) {
      // controlled-R_m gate between qubit j (target) and qubit k (control),
      // angle 2 * \pi / 2^(k - j + 1)
      int m = k - j + 1;
      double phi = 2.0 * M_PI / (double)(1LL << m);

      qstate_apply_controlled_phase(psi, n_qubits, qubits[k], qubits[j], phi);
    }
  }

  for (int j = 0; j < n_q / 2; j++) {
    qstate_apply_swap(psi, n_qubits, qubits[j], qubits[n_q - 1 - j]);
  }
}

void qft_apply_inverse(cvector_t *psi, int n_qubits, const int *qubits,
                       int n_q) {
  if (!psi || !qubits || n_q < 1) {
    return;
  }

  for (int j = 0; j < n_q / 2; j++) {
    qstate_apply_swap(psi, n_qubits, qubits[j], qubits[n_q - 1 - j]);
  }

  for (int j = n_q - 1; j >= 0; j--) {
    for (int k = n_q - 1; k > j; k--) {
      int m = k - j + 1;
      double phi = -2.0 * M_PI / (double)(1LL << m);

      qstate_apply_controlled_phase(psi, n_qubits, qubits[k], qubits[j], phi);
    }

    qstate_apply_gate1(psi, n_qubits, qubits[j], hadamard_gate);
  }
}
