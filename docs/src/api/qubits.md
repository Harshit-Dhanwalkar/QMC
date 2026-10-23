# `qubits`

> Generated from `physics/qubits.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Minimal multi-qubit state-vector simulator: single-qubit gates and CNOT.

**Category:** Quantum computing & information

```c
#include <qmc/physics/qubits.h>
```

**Start with:** [`qstate_alloc`](#qstate_alloc), [`qstate_probability`](#qstate_probability), [`qstate_apply_gate1`](#qstate_apply_gate1), [`qstate_apply_controlled_u`](#qstate_apply_controlled_u), [`qstate_apply_cnot`](#qstate_apply_cnot)

- Guide: [physics/qubits.md](../physics/qubits.md)
- Example: `examples/eg_18_qubits.c`
- Test (every call with expected results): `tests/test_qubits.c`
- Runtime help: `qmc_help(stdout, "qubits");` or `qmc-help qubits`

## Overview

```text
NOTE: Minimal multi-qubit state vector substrate: single-qubit gates + one
entangling two-qubit gate (CNOT) for universal quantum computation
(Solovay-Kitaev: arbitrary single-qubit rotations + one entangler generate
any N-qubit unitary to arbitrary precision).
```

## Data

```c
extern const complex_t hadamard_gate[4];
```

## Functions (15)

[`qstate_alloc`](#qstate_alloc), [`qstate_probability`](#qstate_probability), [`qstate_apply_gate1`](#qstate_apply_gate1), [`qstate_apply_controlled_u`](#qstate_apply_controlled_u), [`qstate_apply_cnot`](#qstate_apply_cnot), [`qstate_apply_swap`](#qstate_apply_swap), [`qstate_apply_controlled_phase`](#qstate_apply_controlled_phase), [`qstate_apply_controlled_unitary`](#qstate_apply_controlled_unitary), [`rx_gate`](#rx_gate), [`ry_gate`](#ry_gate), [`rz_gate`](#rz_gate), [`qstate_reduced_density_single`](#qstate_reduced_density_single), [`von_neumann_entropy_2x2`](#von_neumann_entropy_2x2), [`qstate_measure`](#qstate_measure), [`qstate_measure_qubit`](#qstate_measure_qubit)

### `qstate_alloc`

```c
cvector_t *qstate_alloc(int n_qubits);
```

```text
Allocate an n-qubit state initialized to |00...0>, n >= 1 and small enough
that 2^n doubles fit in memory.
```

### `qstate_probability`

```c
double qstate_probability(const cvector_t *psi, int index);
```

```text
P(measuring computational basis state `index`) = |amplitude|^2
```

### `qstate_apply_gate1`

```c
void qstate_apply_gate1(cvector_t *psi, int n_qubits, int target,
                        const complex_t gate[4]);
```

```text
Apply an arbitrary single-qubit gate to qubit `target` (0-indexed,
0=leftmost) of an n_qubits-qubit state \psi (length 2^n_qubits).
```

### `qstate_apply_controlled_u`

```c
void qstate_apply_controlled_u(cvector_t *psi, int n_qubits, int control,
                               int target, const complex_t U[4]);
```

```text
NOTE: CNOT is special case U = \sigma_x (from angular.c)
```

### `qstate_apply_cnot`

```c
void qstate_apply_cnot(cvector_t *psi, int n_qubits, int control, int target);
```

```text
CNOT(control, target) = qstate_apply_controlled_u with U = \sigma_x
```

### `qstate_apply_swap`

```c
void qstate_apply_swap(cvector_t *psi, int n_qubits, int qa, int qb);
```

```text
SWAP(qa, qb): exchange the states of two qubits.
```

### `qstate_apply_controlled_phase`

```c
void qstate_apply_controlled_phase(cvector_t *psi, int n_qubits, int control,
                                   int target, double phi);
```

```text
Controlled-phase: multiply amplitude by \exp^{i * \phi} whenever both
`control` and `target` are |1>. This is CPHASE/CR_k, the building block of
the QFT's controlled-rotation ladder (a diagonal 2-qubit gate, so it commutes
with itself).
```

### `qstate_apply_controlled_unitary`

```c
void qstate_apply_controlled_unitary(cvector_t *psi, int n_qubits, int control,
                                     const int *targets, int n_targets,
                                     const cmatrix_t *U);
```

```text
General controlled-U primitive for an m-qubit (not just single-qubit)
unitary: apply the 2^m x 2^m unitary matrix U to the `n_targets` qubits
listed in `targets` (in the same MSB-first index convention as every other
qstate_* function; `targets[0]` is the most-significant of the m-qubit
sub-register U acts on), but only on the branch where `control` is |1>.
`n_targets` up to ~5 in practice (block = 2^n_targets amplitudes gathered per
basis-state group, stack-allocated). Needed for Quantum Phase Estimation's
controlled-U^(2^k) step when the eigenstate register spans more than one
qubit (qstate_apply_controlled_u only covers the single-target-qubit case).
```

### `rx_gate`

```c
void rx_gate(double theta, complex_t out[4]);
```

```text
Parameterized single-qubit rotation gates, writing 2x2 gate matrix
(row-major, same [4]-array convention as qstate_apply_gate1's `gate`
argument) into `out`:
RX(\theta) = [[\cos(\theta / 2), -i * \sin(\theta / 2)],
              [-i * \sin(\theta / 2), \cos(\theta / 2)]]
RY(\theta) = [[\cos(\theta / 2), -\sin(\theta / 2],
              [\sin(\theta / 2), \cos(\theta / 2)]
RZ(\theta) = [[\exp(-i * \theta/2), 0],
              [0, \exp(i * \theta / 2)]]

Generators \exp(-i * \theta/2 * \sigma_{x,y,z}), each has period 4*\pi in
\theta. Used as parameterized building blocks for variational ansatz
circuits.
```

### `ry_gate`

```c
void ry_gate(double theta, complex_t out[4]);
```

_No description in the header yet._

### `rz_gate`

```c
void rz_gate(double theta, complex_t out[4]);
```

_No description in the header yet._

### `qstate_reduced_density_single`

```c
cmatrix_t *qstate_reduced_density_single(const cvector_t *psi, int n_qubits,
                                         int qubit);
```

```text
Reduced density matrix of single qubit, tracing out all others:
  \rho_{ab} = \sum_{rest} conj(\psi[rest,qubit=a]) * \psi[rest,qubit=b]

Returns a 2x2 cmatrix_t.
```

### `von_neumann_entropy_2x2`

```c
double von_neumann_entropy_2x2(cmatrix_t *rho);
```

```text
Von Neumann entropy :
  S = -\sum_i \lambda_i \log2(\lambda_i)
of a 2x2 Hermitian density matrix, computed via closed-form eigenvalues of
2x2 Hermitian matrix (\lambda = (tr \pm \sqrt(tr^2 - 4 * det)) / 2)

NOTE: Not solved via general eigensolver, so entropy has no dependency on
eigen_t's internal.

Returns entropy in bits (log base 2): 0 for pure/unentangled qubit, 1 for
maximally mixed (maximally entangled) qubit.
```

### `qstate_measure`

```c
int qstate_measure(cvector_t *psi, double u);
```

```text
Projective measurement in computational basis, on pure state vector.
Samples an outcome index in [0, \psi->n) with probability |\psi[outcome]|^2,
using caller-supplied uniform random number u in [0,1) so RNG/seeding stays
under caller control. Collapses \psi in place to pure basis state |outcome>
(\psi[outcome] = 1, every other amplitude = 0).

Returns sampled outcome index, or -1 on invalid input.
```

### `qstate_measure_qubit`

```c
int qstate_measure_qubit(cvector_t *psi, int n_qubits, int target, double u);
```

```text
Projective measurement of single qubit `target` (0-indexed, 0=leftmost) in
computational basis, leaving other n_qubits-1 qubits in renormalized
superposition survives outcome

P(outcome=0) = sum over all basis states i with bit `target` = 0 of
|\psi[i]|^2. Samples outcome in {0,1} using caller-supplied uniform random
`u` in [0,1). Collapses \psi in place: amplitudes inconsistent with sampled
outcome are zeroed, survivors are rescaled by 1/\sqrt(P(outcome)) so state
remains normalized.

\psi->n must equal 2^n_qubits (caller-supplied n_qubits, since state vector
itself carries no qubit-count metadata).

Returns sampled outcome (0 or 1), or -1 on invalid input (NULL \psi, target
out of [0, n_qubits), or psi->n != 2^n_qubits).
```
