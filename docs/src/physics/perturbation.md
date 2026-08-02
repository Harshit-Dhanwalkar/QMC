# Perturbation Theory

Perturbation theory is a method for finding approximate solutions to quantum systems that cannot be solved exactly. The Hamiltonian is split into an exactly solvable part plus a small perturbation.

## Non-Degenerate Perturbation Theory

The Hamiltonian: $H = H_0 + \lambda V$

The unperturbed eigenstates: $H_0|n\rangle = E_n^{(0)}|n\rangle$

### First Order Correction

$$
E_n^{(1)} = \langle n|V|n\rangle
$$

### Second Order Correction

$$
E_n^{(2)} = \sum_{k \neq n} \frac{|\langle k|V|n\rangle|^2}{E_n^{(0)} - E_k^{(0)}}
$$

## Implementation

`physics/perturbation.h`:

```c
typedef struct {
  double E0; // unperturbed energy
  double E1; // first-order correction
  double E2; // second-order correction
} perturb_result_t;

perturb_result_t perturb_nondeg(const eigen_t *unperturbed, int state_index,
                                const cmatrix_t *V_pert, double tol);
```

`tol` guards the second-order sum against near-degenerate denominators $E_n^{(0)} - E_k^{(0)} \approx 0$ - presumably skipping or otherwise handling terms below that threshold, though the exact behavior needs `perturbation.c` to confirm.

```c
eigen_t *unperturbed = /* solved H0 eigenstates, e.g. via tridiag_eigh */;
cmatrix_t *V_pert = /* build <i|V|j> matrix in the unperturbed basis */;

perturb_result_t result = perturb_nondeg(unperturbed, n, V_pert, 1e-10);
double E_total = result.E0 + lambda * result.E1 + lambda * lambda * result.E2;
```

## Example: Anharmonic Oscillator

For the anharmonic oscillator with perturbation $V = \lambda x^4$, `V_pert` is built as the matrix of $\langle i|x^4|j\rangle$ in the harmonic-oscillator eigenbasis (obtained from `tridiag_eigh`, see [Linear Algebra Core](../internals/linalg.md)), then passed to
`perturb_nondeg` as above.

### Degenerate Perturbation Theory

When two or more unperturbed states share an energy, non-degenerate perturbation theory's second-order sum blows up ($E_n^{(0)} - E_k^{(0)} = 0$). The fix is to diagonalize the perturbation within the degenerate subspace first:

```c
eigen_t *perturb_degenerate(const double *energies, const cmatrix_t *V_pert,
                            const int *degeneracy_indices, int deg_size);
```

`degeneracy_indices` lists which indices (into `energies` and `V_pert`) belong to the degenerate subspace; the returned `eigen_t` has size `deg_size` and gives the correct zeroth-order states and first-order energy splittings within that subspace.

### Fermi's Golden Rule

For time-dependent perturbations, the transition rate from state $i$ to $f$ is:

$$
W_{i \rightarrow f}=\frac{2\pi}{\hbar}\vert{}\langle f\vert{}V\vert{}i \rangle\vert{}^2 \rho(E_f)
$$

```c
double fermi_golden_rate(const cmatrix_t *V_pert, int i, int f, double rho_E);
```

`V_pert` is the perturbation operator in the relevant basis; `rho_E` is the density of final states at $E_f \approx E_i$, computed separately and passed in rather than derived internally.
