# `lindblad`

> Generated from `physics/lindblad.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Lindblad (GKSL) master equation for open quantum systems.

**Category:** Two-level, open systems & quantum optics

```c
#include <qmc/physics/lindblad.h>
```

**Start with:** [`density_from_pure_state`](#density_from_pure_state), [`density_purity`](#density_purity), [`density_von_neumann_entropy`](#density_von_neumann_entropy), [`lindblad_rhs`](#lindblad_rhs), [`lindblad_step_rk4`](#lindblad_step_rk4)

- Guide: [physics/lindblad.md](../physics/lindblad.md)
- Example: `examples/eg_19_lindblad.c`
- Test (every call with expected results): `tests/test_lindblad.c`
- Runtime help: `qmc_help(stdout, "lindblad");` or `qmc-help lindblad`

## Overview

```text
Lindblad (GKSL) equation for open quantum systems:

  d\rho/dt = -i [H, \rho] + \sum_k D[L_k](\rho)
  D[L](\rho) = L \rho L^dag - 1/2 { L^dag L, \rho }

Natural units (\hbar = m = 1)
\rho is an NxN Hermitian density matrix (cmatrix_t,
Tr(\rho) = 1); H is NxN Hermitian; each L_k is an NxN jump operator
```

## Functions (11)

[`density_from_pure_state`](#density_from_pure_state), [`density_purity`](#density_purity), [`density_von_neumann_entropy`](#density_von_neumann_entropy), [`lindblad_rhs`](#lindblad_rhs), [`lindblad_step_rk4`](#lindblad_step_rk4), [`lindblad_evolve`](#lindblad_evolve), [`embed_single_qubit_op`](#embed_single_qubit_op), [`lindblad_amplitude_damping_op`](#lindblad_amplitude_damping_op), [`lindblad_dephasing_op`](#lindblad_dephasing_op), [`lindblad_bitflip_op`](#lindblad_bitflip_op), [`density_measure_computational_basis`](#density_measure_computational_basis)

### `density_from_pure_state`

```c
cmatrix_t *density_from_pure_state(const cvector_t *psi);
```

```text
rho = |\psi><\psi| (outer product), an NxN density matrix from an
N-dimensional pure state vector
```

### `density_purity`

```c
double density_purity(const cmatrix_t *rho);
```

```text
Purity Tr(\rho^2) = \sum_{ij} |\rho_{ij}|^2
1 for pure state, < 1 for a mixed state, 1/N for maximally mixed state
```

### `density_von_neumann_entropy`

```c
double density_von_neumann_entropy(const cmatrix_t *rho);
```

```text
Von Neumann entropy S = -Tr(\rho log2 \rho) of general NxN Hermitian density
matrix

Returns entropy in bits: 0 for a pure state, log2(N) for maximally mixed
```

### `lindblad_rhs`

```c
cmatrix_t *lindblad_rhs(const cmatrix_t *H, const cmatrix_t *rho, cmatrix_t **L,
                        int n_ops);
```

```text
Right-hand side of Lindblad equation at current \rho:
  -i[H,\rho] + \sum_k(L_k \rho L_k^\dagger - 1/2 {L_k^\dagger L_k, \rho})

L is an array of n_ops jump operators (each NxN, N = \rho->nrows)
Pass n_ops = 0 (L may be NULL) for closed-system (unitary) evolution

Returns allocated NxN cmatrix_t (caller frees), or NULL on dimension mismatch
/ allocation failure
```

### `lindblad_step_rk4`

```c
int lindblad_step_rk4(const cmatrix_t *rho, const cmatrix_t *H, cmatrix_t **L,
                      int n_ops, double dt);
```

```text
Advance \rho by one step of size dt using 4th-order Runge-Kutta applied to
lindblad_rhs

Returns 0 on success, -1 on invalid input (dimension mismatch, etc)
```

### `lindblad_evolve`

```c
int lindblad_evolve(const cmatrix_t *rho, const cmatrix_t *H, cmatrix_t **L,
                    int n_ops, double dt, int steps);
```

```text
Evolve \rho for `steps` of size dt via repeated lindblad_step_rk4

Returns 0 on success, -1 on first failed step (\rho may be partially evolved
up to that point)
```

### `embed_single_qubit_op`

```c
cmatrix_t *embed_single_qubit_op(const complex_t op[4], int n_qubits,
                                 int target);
```

```text
Embed a single-qubit 2x2 operator `op` (row-major: op[0]=<0|op|0>,
op[1]=<0|op|1>, op[2]=<1|op|0>, op[3]=<1|op|1>) acting on qubit `target`
of an n_qubits-qubit system into full 2^n_qubits x 2^n_qubits Hilbert
space (op on `target`, identity on every other qubit)

Returns a newly allocated (2^n_qubits)x(2^n_qubits) cmatrix_t, or NULL on
invalid input
```

### `lindblad_amplitude_damping_op`

```c
cmatrix_t *lindblad_amplitude_damping_op(int n_qubits, int target,
                                         double gamma);
```

```text
Common single-qubit noise channels, embedded into an n_qubits-qubit Hilbert
space and acting on qubit `target`, as an L_k. \gamma is channel's rate
(1 / time, natural units)

 - amplitude damping (T1 decay): L = \sqrt(\gamma) * \sigma_minus,
   \sigma_minus = |0><1| (population decays from |1> to |0>)
 - pure dephasing (T2):          L = \sqrt(\gamma/2) * \sigma_z
 - bit flip:                     L = \sqrt(\gamma) * \sigma_x
```

### `lindblad_dephasing_op`

```c
cmatrix_t *lindblad_dephasing_op(int n_qubits, int target, double gamma);
```

_No description in the header yet._

### `lindblad_bitflip_op`

```c
cmatrix_t *lindblad_bitflip_op(int n_qubits, int target, double gamma);
```

_No description in the header yet._

### `density_measure_computational_basis`

```c
int density_measure_computational_basis(cmatrix_t *rho, double u);
```

```text
Projective measurement in computational basis:
Samples an outcome index in [0, \rho->nrows) according to diagonal
populations \rho[i][i].re, using caller-supplied uniform random number `u` in
[0,1) (so RNG/seeding stays under caller control) Collapses \rho in place to
pure-state projector |outcome><outcome| for sampled outcome

Returns sampled outcome index, or -1 on invalid input (\rho not square)
```
