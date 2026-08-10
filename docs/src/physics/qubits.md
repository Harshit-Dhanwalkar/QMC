# Qubits & Quantum Computation

Qubit operations represent single- and multi-qubit transformations using density matrices and state vectors.

A minimal multi-qubit state-vector substrate: single-qubit gates plus one entangling two-qubit gate (CNOT). By the Solovay-Kitaev theorem, arbitrary single-qubit rotations plus one entangler are sufficient to generate any $N$-qubit unitary to arbitrary precision - so this small set is computationally universal. Implemented in `physics/qubits.h`.

## Single Qubit State

$$
|\psi\rangle = \alpha |0\rangle + \beta |1\rangle = \cos\left(\frac{\theta}{2}\right)|0\rangle + e^{i\phi}\sin\left(\frac{\theta}{2}\right)|1\rangle
$$

## Standard Quantum Gates

- **Pauli Gates**:
  $$
  X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad Y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad Z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}
  $$
- **Hadamard Gate**:
  $$
  H = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}
  $$
- **CNOT Gate**:
  $$
  \text{CNOT} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \end{pmatrix}
  $$

## State representation

An $n$-qubit state is a `cvector_t` of length $2^n$ - no separate tensor or register type:

```c
// |00...0>, n >= 1 and small enough that 2^n doubles fit in memory
cvector_t *qstate_alloc(int n_qubits);

// P(measuring computational basis state `index`) = |amplitude|^2
double qstate_probability(const cvector_t *psi, int index);
```

## Gates

```c
void qstate_apply_gate1(cvector_t *psi, int n_qubits, int target,
                        const complex_t gate[4]);
```

Applies an arbitrary single-qubit gate (as a flat row-major 2x2, same convention as the Pauli matrices in [Angular Momentum & Spin](angular_momentum.md)) to qubit `target` (0-indexed, 0 = leftmost) of an `n_qubits`-qubit state.

```c
extern const complex_t hadamard_gate[4]; // (1/\sqrt2) [[1,1],[1,-1]]
```

### Controlled gates

```c
void qstate_apply_controlled_u(cvector_t *psi, int n_qubits, int control,
                               int target, const complex_t U[4]);

void qstate_apply_cnot(cvector_t *psi, int n_qubits, int control, int target);
```

`qstate_apply_controlled_u` applies `U` to `target` only when `control` is $|1\rangle$ (`control != target` required). CNOT is the special case $U = \sigma_x$ (from `angular.c`) - `qstate_apply_cnot` is just a thin convenience wrapper around `qstate_apply_controlled_u`.

```c
cvector_t *psi = qstate_alloc(2);             // |00>
qstate_apply_gate1(psi, 2, 0, hadamard_gate); // qubit 0 -> superposition
qstate_apply_cnot(psi, 2, 0, 1);              // entangle: Bell state
```

## Entanglement

```c
// \rho_ab = \sum_{rest} conj(\psi[rest,qubit=a]) * \psi[rest,qubit=b]
cmatrix_t *qstate_reduced_density_single(const cvector_t *psi, int n_qubits,
                                         int qubit);

/* Von Neumann entropy in bits: 0 for pure/unentangled, 1 for maximally
   mixed (maximally entangled) qubit */
double von_neumann_entropy_2x2(cmatrix_t *rho);
```

`von_neumann_entropy_2x2` uses the closed-form eigenvalues of a 2x2 Hermitian matrix ($\lambda = (\text{tr} \pm \sqrt{\text{tr}^2 - 4\det})/2$) rather than a general eigensolver - the header explicitly notes this has no dependency on `eigen_t`'s internals, so it's independent of any changes to [the eigensolvers](../internals/linalg.md).

```c
cmatrix_t *rho0 = qstate_reduced_density_single(psi, 2, 0);
double S = von_neumann_entropy_2x2(rho0); /* 1.0 for the Bell state above */
```

## Measurement

```c
/* Samples an outcome in [0, psi->n) with probability |psi[outcome]|^2, using
   caller-supplied uniform u in [0,1) so RNG/seeding stays under caller
   control. Collapses psi in place. Returns outcome index, or -1 on invalid
   input. */
int qstate_measure(cvector_t *psi, double u);
```

Deterministic given `u` - the caller owns the RNG, which keeps measurement reproducible for tests.

## Running the Example

```sh
./build/eg_18_qubits
```
