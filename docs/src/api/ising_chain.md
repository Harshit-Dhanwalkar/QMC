# `ising_chain`

> Generated from `physics/ising_chain.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Transverse-field Ising chain: exact diagonalisation on a ring and the exact Jordan-Wigner solution.

**Category:** Condensed matter & many-body

```c
#include <qmc/physics/ising_chain.h>
```

**Start with:** [`ising_hamiltonian`](#ising_hamiltonian), [`ising_z2_hamiltonian`](#ising_z2_hamiltonian), [`ising_exact_ground_energy_per_site`](#ising_exact_ground_energy_per_site), [`ising_reduced_density_matrix`](#ising_reduced_density_matrix), [`ising_entanglement_entropy`](#ising_entanglement_entropy)

- Example: `examples/eg_63_ising_chain.c`
- Test (every call with expected results): `tests/test_ising_chain.c`
- Runtime help: `qmc_help(stdout, "ising_chain");` or `qmc-help ising_chain`

## Overview

```text
Exact diagonalization for spin-1/2 transverse-field Ising model (TFIM) on a
ring, plus its exact Jordan-Wigner + Bogoliubov ground-state energy

H = -J \sum_j \sigma^z_j \sigma^z_{j+1} - h \sum_j \sigma^x_j

NOTE: Convention: site states are packed into bits of a uint32_t (bit j = 0
-> spin up, \sigma^z_j = +1; bit j = 1 -> spin down, \sigma^z_j = -1). This
limits N <= 30, far beyond what an explicit 2^N x 2^N (or even 2^(N-1) x
2^(N-1)) sparse matrix can hold anyway
```

## Functions (7)

[`ising_hamiltonian`](#ising_hamiltonian), [`ising_z2_hamiltonian`](#ising_z2_hamiltonian), [`ising_exact_ground_energy_per_site`](#ising_exact_ground_energy_per_site), [`ising_reduced_density_matrix`](#ising_reduced_density_matrix), [`ising_entanglement_entropy`](#ising_entanglement_entropy), [`ising_hamiltonian_dense`](#ising_hamiltonian_dense), [`ising_loschmidt_echo`](#ising_loschmidt_echo)

### `ising_hamiltonian`

```c
sparse_matrix_t *ising_hamiltonian(int N, double J, double h, int pbc);
```

```text
Full (unreduced) TFIM Hamiltonian over all 2^N computational basis states, in
CSR form (reuses core/sparse.h so it plugs directly into lanczos_eigs /
lanczos_tridiagonalize). pbc selects periodic (ring) vs open boundary
conditions

Returns NULL for N <= 0, N > 30, or allocation failure
```

### `ising_z2_hamiltonian`

```c
sparse_matrix_t *ising_z2_hamiltonian(int N, double J, double h, int pbc,
                                      int parity);
```

```text
Z2 (global spin-flip) parity-reduced Hamiltonian, dimension 2^(N-1). parity
must be exactly +1 or -1. Row/column a corresponds to a-th smallest
computational basis state s satisfying s < (s XOR mask) (i.e. canonical,
"smaller" representative of each flip-pair), in ascending order

TFIM ground state is always in +1 sector for h,J >= 0, but both sectors are
exposed since gap between them controls e.g. (exponentially small, for finite
N in ordered phase) tunneling splitting between 2 symmetry-broken
ferromagnetic states

Returns NULL for N <= 0, N > 30, parity not in {+1,-1}, or allocation failure
```

### `ising_exact_ground_energy_per_site`

```c
double ising_exact_ground_energy_per_site(int N, double J, double h);
```

```text
Exact ground-state energy per site of periodic transverse-field Ising model,
via Jordan-Wigner fermionization into free fermions followed by a Bogoliubov
transformation:

  E0/N = -(1/N) \sum_k \epsilon_k / ...

Where
 \epsilon_k = 2 * \sqrt(J^2 + h^2 - 2 * J * h * \cos(k)) are Bogoliubov
quasiparticle energies. This is completely independent of ED machinery above
provided purely as a validation reference

Valid for any N >= 1; J, h may be any real values (h=0 or J=0 are physically
fine limits, exactly solvable and non-degenerate whenever h != J)
```

### `ising_reduced_density_matrix`

```c
cmatrix_t *ising_reduced_density_matrix(const cmatrix_t *psi, int N, int L_A);
```

```text
Reduced density matrix of "A" subsystem consisting of sites 0..L_A-1 (low L_A
bits of state index), obtained by tracing out remaining N-L_A sites ("B")
from a pure state `\psi` (a computational-basis state vector of dimension
2^N, e.g. a ground-state eigenvector from lanczos_eigs on ising_hamiltonian's
output)

  \rho_A[a][a'] = \sum_b \psi[a | (b<<L_A)] * conj(\psi[a' | (b<<L_A)])

Returns a dense 2^L_A x 2^L_A Hermitian matrix (Tr(\rho_A) = 1 whenever \psi
is normalized). L_A=0 or L_A=N both correctly give a trivial 1x1 matrix [1]
(respectively: subsystem A empty, or A is whole system - pure global state
has zero entanglement with "nothing"). Returns NULL for N <= 0, N > 30, L_A <
0, L_A > N, or allocation failure
```

### `ising_entanglement_entropy`

```c
double ising_entanglement_entropy(const cmatrix_t *psi, int N, int L_A);
```

```text
Von Neumann entanglement entropy S = -Tr(rho_A log rho_A) of subsystem
consisting of sites 0..L_A-1, for a pure state `\psi`. Diagonalizes
ising_reduced_density_matrix(psi, N, L_A) internally (dense complex Hermitian
eigensolver) and sums -p * \ln(p) over its eigenvalues, skipping any at or
below a small numerical floor (exact zero eigenvalues would otherwise make
log(0) blow up, and floating-point noise can occasionally push a true zero
slightly negative)

Two exact, parameter-independent identities this satisfies for any valid pure
state : S(L_A) = S(N-L_A) (a pure global state's subsystem and complement
always have equal entanglement entropy), and S(0)=S(N)=0 (trivial subsystem
has nothing to be entangled with)

Returns NAN for same invalid inputs ising_reduced_density_matrix rejects, or
if diagonalization fails
```

### `ising_hamiltonian_dense`

```c
cmatrix_t *ising_hamiltonian_dense(int N, double J, double h, int pbc);
```

```text
Dense N x N counterpart of ising_hamiltonian, for small-N calculations that
need a FULL eigendecomposition (e.g. quench dynamics below) rather than just
lowest few eigenpairs Lanczos gives. Same Hamiltonian, same conventions;
returns NULL for same invalid inputs. Intended for N small enough that
dim=2^N is comfortably diagonalizable dense
```

### `ising_loschmidt_echo`

```c
double ising_loschmidt_echo(const eigen_t *eig_f, const cmatrix_t *psi_i,
                            double t);
```

```text
Loschmidt echo L(t) = |<\psi_i| \exp(-i * H_f * t) |\psi_i>|^2 following a
sudden ("quantum quench") change of transverse field, given H_f's full
eigendecomposition eig_f (e.g. from cmatrix_eigh_complex on
ising_hamiltonian_dense's output) and pre-quench state \psi_i (typically
ground-state eigenvector of a different hamiltonian, at pre-quench field
value h_i != h_f)

Exact identities that hold for any valid Hermitian H_f and pure state \psi_i
not specific to TFIM : L(0) = 1 always; if \psi_i happens to already be an
eigenvector of H_f (trivial "no quench" case), L(t) = 1 for all t; and L(t)
admits short-time expansion

  L(t) = 1 - t^2 * Var(H_f)_{\psi_i} + O(t^4)

Where
 Var(H_f)_{\psi_i} = <\psi_i|H_f^2|\psi_i> - <\psi_i|H_f|\psi_i>^2

Returns NAN if eig_f or \psi_i is NULL, or if their dimensions don't match
```
