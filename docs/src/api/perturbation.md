# `perturbation`

> Generated from `physics/perturbation.h` by `tools/gen_api.py`. To change this page, edit the header's comments or `tools/api_meta.py`.

Time-independent perturbation theory (non-degenerate and degenerate) and Fermi's golden rule.

**Category:** Atomic structure & perturbation theory

```c
#include <qmc/physics/perturbation.h>
```

**Start with:** [`perturb_nondeg`](#perturb_nondeg), [`perturb_degenerate`](#perturb_degenerate), [`fermi_golden_rate`](#fermi_golden_rate)

- Guide: [physics/perturbation.md](../physics/perturbation.md)
- Example: `examples/eg_10_perturbation.c`
- Test (every call with expected results): `tests/test_perturbation.c`
- Runtime help: `qmc_help(stdout, "perturbation");` or `qmc-help perturbation`

## Types

### `perturb_result_t`

```text
Non-degenerate perturbation theory (first and second order)
```

```c
typedef struct {
  double E0; // unperturbed energy
  double E1; // first-order correction
  double E2; // second-order correction

  /* HACK: diagnostics for E2: perturb_nondeg's second-order sum omits any term
   * whose energy denominator |E_i - E_k| falls below `tol` (a near-degenerate
   * state k), since non-degenerate PT's 1/(E_i-E_k) formula diverges there
   *
   * Where
   *   n_terms_skipped_near_degenerate: how many k were omitted (0 for a clean
   *                                    non-degenerate calculation)
   *   min_skipped_denominator        : min(|E_i-E_k|) among omitted terms (only
   *                                    meaningful when
   *                                    n_terms_skipped_near_degenerate > 0),
   *                                    closer to 0, more E2 likely
   * underestimates true correction.
   *
   * When n_terms_skipped_near_degenerate > 0, use perturb_degenerate on
   * near-degenerate subspace instead of trusting E2 from this call */
  int n_terms_skipped_near_degenerate;
  double min_skipped_denominator;
} perturb_result_t;
```

## Functions (3)

[`perturb_nondeg`](#perturb_nondeg), [`perturb_degenerate`](#perturb_degenerate), [`fermi_golden_rate`](#fermi_golden_rate)

### `perturb_nondeg`

```c
perturb_result_t perturb_nondeg(const eigen_t *unperturbed, int state_index,
                                const cmatrix_t *V_pert, double tol);
```

_No description in the header yet._

### `perturb_degenerate`

```c
eigen_t *perturb_degenerate(const double *energies, const cmatrix_t *V_pert,
                            const int *degeneracy_indices, int deg_size);
```

```text
Degenerate perturbation theory: diagonalize V_pert in degenerate subspace

Returns new eigen_t with energies (eigenvalues) and eigenvectors (in that
subspace)
Subspace is defined by indices of states with same (or close) energy

Where
  energies          : array of unperturbed energies (size n)
  V_pert            : perturbation matrix (size n x n)
  degeneracy_indices: array of indices that belong to degenerate subspace

Returns eigen_t of effective Hamiltonian (size = subspace size)
```

### `fermi_golden_rate`

```c
double fermi_golden_rate(const cmatrix_t *V_pert, int i, int f, double rho_E);
```

```text
Time-dependent perturbation: Fermi's golden rule
transition rate from initial state |i> to final states with density of states
\rho(E)

Where
  V_pert  : perturbation operator (matrix element <f|V|i>)
  \rho_E  : density of states at energy E_f ~ E_i
```
