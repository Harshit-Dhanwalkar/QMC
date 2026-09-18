# Complete Active Space Self-Consistent Field (CASSCF)

Unlike CISD/FCI, which diagonalize a Hamiltonian built from a _fixed_ set of
orbitals, CASSCF simultaneously optimizes the CI coefficients _and_ the
molecular orbitals themselves. Implemented in `physics/casscf.h`/
`physics/casscf.c`, on top of `molecular_hf.c` and `molecular_integrals.c`.

## Active space partition

The $n_{\text{basis}}$ spatial orbitals are split into three groups, the
same partition every quantum chemistry package uses for CASSCF:

- **frozen** ($n_{\text{frozen}}$): always doubly occupied, excluded from
  the CI expansion.
- **active** ($n_{\text{active}}$): a full CI is run within this space.
- **virtual** (the remainder, $n_{\text{basis}} - n_{\text{frozen}} -
  n_{\text{active}}$): always empty.

## Algorithm

1. **Active-space CI**: build the frozen-core-effective Hamiltonian over
   just the $n_{\text{frozen}}+n_{\text{active}}$ orbitals, diagonalize in
   the $n_{\text{electrons\_active}}$ sector, and construct the 1- and
   2-particle reduced density matrices (RDMs) from the resulting CI vector
   via direct bit-manipulation fermionic operators.
2. **Extend to the full orbital space**: core orbitals contribute a fixed
   closed-shell density block, active contributes the CI density, virtual
   is zero (standard MCSCF density partition).
3. **Orbital gradient**: rather than an analytically-derived generalized
   Fock matrix, $\partial E/\partial\kappa_{pq}$ for every independent
   orbital-rotation generator $\kappa_{pq}$ is computed by central finite
   difference, holding the RDMs from step 1 fixed and only recomputing the
   energy expression
   $E = E_{\text{nuc}} + \sum h_{pq} D_{pq} + \frac12\sum (pq|rs)\, d_{pqrs}$
   under a small orbital rotation (no re-diagonalization needed). The
   Hessian diagonal is estimated the same way (3-point finite difference),
   giving an approximate Newton step per generator; the denominator is
   floored in absolute value (not sign-preserved) to guarantee a valid
   descent direction even where the true curvature is tiny or slightly
   negative.
4. **Line search**: backtrack on the actual (re-diagonalized) active-space
   CI energy, then rotate $C \leftarrow C\, e^{\kappa}$ using a real matrix
   exponential (scaling-and-squaring + Taylor series).
5. Repeat from step 1 with the rotated orbitals until the gradient norm
   falls below `conv_tol`.

## Implementation

```c
typedef struct {
  double total_energy; // electronic + nuclear repulsion
  int n_spatial;
  int n_frozen;
  int n_active;
  int n_electrons_active;
  cmatrix_t *C;     // optimized n_spatial x n_spatial MO coefficients
  int iterations;   // orbital-rotation steps taken
  int converged;    // 1 if final gradient norm < conv_tol, else 0
  double grad_norm; // final orbital-gradient norm
} casscf_result_t;

casscf_result_t *casscf_run(basis_function_t **basis, int n_basis,
                            const molecule_t *mol, const cmatrix_t *C_initial,
                            int n_frozen, int n_active, int n_electrons_active,
                            double conv_tol, int max_iter);

void casscf_result_free(casscf_result_t *res);
```

`C_initial` is a starting MO coefficient matrix, typically a converged RHF
result's `C`: CASSCF is a local optimization and, like any MCSCF method,
needs a reasonable starting guess. Requires `0 <= n_frozen`,
`n_active >= 1`, `n_frozen + n_active <= n_basis`, and
`0 <= n_electrons_active <= 2*n_active`. Returns `NULL` on invalid input or
allocation failure; always check `converged` before trusting `total_energy`
: a non-converged result still holds the best `C`/energy found within
`max_iter` iterations.

```c
molecular_hf_result_t *rhf = molecular_hf_run(/* ... */);
casscf_result_t *cas = casscf_run(basis, n_basis, &mol, rhf->C,
                                  /*n_frozen=*/1, /*n_active=*/2,
                                  /*n_electrons_active=*/2, 1e-6, 100);
if (cas->converged) {
    printf("CASSCF total energy: %.8f Hartree\n", cas->total_energy);
}
casscf_result_free(cas);
```

## Running the Example

```sh
./build/eg_71_casscf
```
