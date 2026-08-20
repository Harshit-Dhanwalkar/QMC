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
