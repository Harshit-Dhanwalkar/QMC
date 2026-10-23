# `spin_chain`

> Generated from `physics/spin_chain.h` by `tools/gen_api.py`. To change this page, edit header's comments or `tools/api_meta.py`.

Spin-1/2 XXZ chain on a ring: symmetry-adapted exact diagonalisation and Lanczos.

**Category:** Condensed matter & many-body

```c
#include <qmc/physics/spin_chain.h>
```

**Start with:** [`spin_sector_build`](#spin_sector_build), [`spin_sector_find`](#spin_sector_find), [`spin_sector_hamiltonian`](#spin_sector_hamiltonian), [`spin_apply_szq`](#spin_apply_szq), [`spin_dsf_continued_fraction`](#spin_dsf_continued_fraction)

- Example: `examples/eg_62_spin_chain_dsf.c`
- Test (every call with expected results): `tests/test_spin_chain.c`
- Runtime help: `qmc_help(stdout, "spin_chain");` or `qmc-help spin_chain`

## Overview

```text
Translation-symmetry-adapted exact diagonalization for the spin-1/2 XXZ
chain on a periodic ring, plus Lanczos continued-fraction evaluation of
the dynamical structure factor S^{zz}(q, \omega).

H = \sum_j [ Jz S^z_j S^z_{j+1}
           + (Jxy/2)(S^+_j S^-_{j+1}
           + S^-_j S^+_{j+1})
           ]

Convention: site states are packed into the bits of a uint64_t (bit j = 0
-> spin up, S^z_j = +1/2; bit j = 1 -> spin down, S^z_j = -1/2). This
limits N <= 63, far beyond what dense/sparse storage can handle anyway.

NOTE: Only two symmetries are exploited here: total S^z (fixed by choosing
number of up spins, nup) and translation (fixed by choosing a momentum index
k in units of 2 * \pi / N). Reflection (spatial parity) and spin-inversion
are also available as further block-diagonalizations on top of these.
```

## Types

### `spin_sector_t`

```c
typedef struct {
  int N;         /* number of sites */
  int nup;       /* number of up spins (fixes the S^z sector) */
  int k;         /* momentum index, 0 <= k < N, in units of 2*pi/N */
  int dim;       /* number of symmetry-adapted basis states in this sector */
  uint64_t *rep; /* size dim, sorted ascending: representative bitstrings */
  int *period;   /* size dim: R_a, the translation period of rep[a] */
} spin_sector_t;
```

## Functions (10)

[`spin_sector_build`](#spin_sector_build), [`spin_sector_free`](#spin_sector_free), [`spin_sector_find`](#spin_sector_find), [`spin_sector_hamiltonian`](#spin_sector_hamiltonian), [`spin_apply_szq`](#spin_apply_szq), [`spin_dsf_continued_fraction`](#spin_dsf_continued_fraction), [`spin_sector_reflection_matrix`](#spin_sector_reflection_matrix), [`spin_sector_parity_project`](#spin_sector_parity_project), [`spin_sector_inversion_matrix`](#spin_sector_inversion_matrix), [`spin_sector_inversion_project`](#spin_sector_inversion_project)

### `spin_sector_build`

```c
spin_sector_t *spin_sector_build(int N, int nup, int k);
```

```text
Build the (nup, k) symmetry sector for an N-site ring. Enumerates all
C(N,nup) basis states, keeps only translation representatives (the
lexicographically/numerically smallest bitstring in each translation
orbit) whose period R_a is compatible with k (i.e. k*R_a = 0 mod N, the
condition for a nonzero-norm momentum eigenstate to exist on that orbit).

Returns NULL on invalid input (N<=0, N>63, nup<0, nup>N, k<0, k>=N) or
allocation failure. Free with spin_sector_free.
```

### `spin_sector_free`

```c
void spin_sector_free(spin_sector_t *sec);
```

_No description in header yet._

### `spin_sector_find`

```c
int spin_sector_find(const spin_sector_t *sec, uint64_t rep);
```

```text
Binary search for a representative bitstring in a sector's sorted rep[]
array. Returns its index, or -1 if not present in this sector.
```

### `spin_sector_hamiltonian`

```c
sparse_matrix_t *spin_sector_hamiltonian(const spin_sector_t *sec, double Jxy,
                                         double Jz, int pbc);
```

```text
Build the matrix-free XXZ Hamiltonian restricted to one symmetry sector,
returned in CSR form (reusing core/sparse.h so it plugs directly into
lanczos_eigs / lanczos_tridiagonalize). pbc selects periodic (ring) vs
open boundary conditions for the j -> j+1 bonds.

Construction is a two-pass, row-parallel (OpenMP) build: pass 1 counts
nonzeros per representative row independently in parallel and prefix-sums
them into row_ptr; pass 2 fills each row's disjoint slice of col_ind/
values in parallel. Because each row's output slice is known in advance
and never touched by another thread, this needs no locks, no atomics, and
no thread-private scratch buffers to reduce afterward.
```

### `spin_apply_szq`

```c
int spin_apply_szq(const spin_sector_t *src, const cvector_t *psi0, int q_index,
                   spin_sector_t **out_target, cvector_t **out_phi0,
                   double *out_I0);
```

```text
Apply S^z_q = (1 / \sqrt(N)) * \sum_j  * \exp(i * q * j) S^z_j to a state
given by its coefficients (psi0, length src->dim) in the src sector's
symmetry-adapted basis, where q = 2 * \pi * q_index / N.

S^z_j is diagonal in the computational (S^z eigenbasis) basis, so this
operator maps src's momentum k to k_target = (src->k + q_index) mod N
while leaving nup (and hence the S^z sector) fixed. The output is
returned unnormalized, together with its squared norm I0 = <phi0|phi0>,
which is itself the q-resolved total spectral weight (equal-time structure
factor at momentum q).

Returns 0 on success, nonzero on invalid input or allocation failure (in
which case no sector/vector is allocated).
```

### `spin_dsf_continued_fraction`

```c
double spin_dsf_continued_fraction(const double *alpha, const double *beta,
                                   int m, double E0, double I0, double omega,
                                   double eta);
```

```text
Evaluate the Lanczos continued-fraction dynamical structure factor
  S(q, \omega) = -(1 / \pi) * Im[ I0 * <f0| (E0 + \omega
                                                + i * \eta
                                                - H)^{-1}
                                       |f0> ]
from tridiagonal coefficients (alpha, beta) already generated by
lanczos_tridiagonalize on normalized excitation vector f0 = phi0/ \sqrt(I0),
plus separately-stored I0 and ground-state energy E0.

NOTE: \eta > 0 is Lorentzian broadening (finite-size chains have
delta-function poles on the real axis, \eta smooths them into a spectrum
which can be plot on a finite frequency mesh).
```

### `spin_sector_reflection_matrix`

```c
cmatrix_t *spin_sector_reflection_matrix(const spin_sector_t *sec);
```

```text
Builds the reflection operator R as a dense sec->dim x sec->dim unitary
matrix in the SAME translation-symmetry-adapted basis (same rep[]/period[]
ordering) as spin_sector_hamiltonian(sec, ...), i.e. the matrix such that
entry (b,a) is <k,rep_b| R |k,rep_a>.

Returns NULL if sec->k is not 0 and not (N even) N/2 : reflection maps
momentum k to -k, so no single-valued reflection quantum number exists on any
other momentum sector (that isn't a bug to fix, it's the actual symmetry
content of the problem). Also returns NULL on allocation failure.
```

### `spin_sector_parity_project`

```c
cmatrix_t *spin_sector_parity_project(const spin_sector_t *sec,
                                      const sparse_matrix_t *H, int parity,
                                      int *dim_out);
```

```text
A sector Hamiltonian H (as returned by spin_sector_hamiltonian on this same
sec) onto the subspace where the reflection operator above has eigenvalue
`parity` (must be exactly +1 or -1), returning a dense Hamiltonian of the
resulting (generally much smaller) dimension in *dim_out. Diagonalizing this
block instead of the full sector reproduces exactly the subset of sec's
spectrum with that reflection parity.

Returns NULL if sec's momentum has no reflection quantum number, if parity
isn't +-1, or on allocation failure. If the requested parity subspace is
empty, returns a valid 0x0 cmatrix_t with *dim_out = 0 (not NULL : that's
reserved for the error cases above).
```

### `spin_sector_inversion_matrix`

```c
cmatrix_t *spin_sector_inversion_matrix(const spin_sector_t *sec);
```

```text
Builds the spin-inversion operator I as a dense sec->dim x sec->dim unitary
matrix in the same translation-symmetry-adapted basis as
spin_sector_hamiltonian(sec, ...).

Returns NULL if sec is not at half filling (N odd, or nup != N/2) - I would
map this sector to a different (N-nup) sector, so no single-valued inversion
quantum number exists there (not a bug: that's the actual symmetry content of
the problem). Also returns NULL on allocation failure.
```

### `spin_sector_inversion_project`

```c
cmatrix_t *spin_sector_inversion_project(const spin_sector_t *sec,
                                         const sparse_matrix_t *H, int parity,
                                         int *dim_out);
```

```text
Sector Hamiltonian H onto subspace where the inversion operator above has
eigenvalue `parity` (must be exactly +1 or -1), returning a dense Hamiltonian
of the resulting (generally much smaller) dimension in  *dim_out.
Diagonalizing this block instead of the full sector reproduces exactly the
subset of sec's spectrum with that inversion parity.

WARN: Unlike spin_sector_parity_project, this works at any momentum k (as
long as sec is at half filling).
Returns NULL if sec isn't at half filling, if parity isn't +-1, or on
allocation failure. If the requested parity subspace is empty, returns a
valid 0x0 cmatrix_t with *dim_out = 0 (not NULL: that's reserved for the
error cases above).
```
