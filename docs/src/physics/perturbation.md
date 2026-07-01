# Perturbation Theory

Perturbation theory is a method for finding approximate solutions to quantum systems that cannot be solved exactly. The Hamiltonian is split into an exactly solvable part plus a small perturbation.

## Non-Degenerate Perturbation Theory

The Hamiltonian: $H = H_0 + \lambda V$

The unperturbed eigenstates: $H_0|n\rangle = E_n^{(0)}|n\rangle$

### First Order Correction

The first-order energy correction:

$$
E_n^{(1)} = \langle n|V|n\rangle
$$

### Second Order Correction

$$
E*n^{(2)} = \sum*{k \neq n} \frac{|\langle k|V|n\rangle|^2}{E_n^{(0)} - E_k^{(0)}}
$$

## Implementation

```c
perturb_result_t perturb_nondeg(const eigen_t *unperturbed, int state_index,
                                const cmatrix_t *V_pert, double tol) {
    perturb_result_t res = {0.0, 0.0, 0.0};

    double E0 = unperturbed->eigenvalues[state_index];
    res.E0 = E0;

    // First order
    complex_t Vii = CMAT(V_pert, state_index, state_index);
    res.E1 = Vii.re;

    // Second order
    double E2 = 0.0;
    for (int k = 0; k < n; k++) {
        if (k == state_index) continue;
        double denom = E0 - unperturbed->eigenvalues[k];
        if (fabs(denom) < tol) continue;
        complex_t Vki = CMAT(V_pert, k, state_index);
        E2 += c_abs2(Vki) / denom;
    }
    res.E2 = E2;
    return res;
}
```

## Example: Anharmonic Oscillator

The anharmonic oscillator with perturbation $V=\lambda x^4$ :

```c
// Build perturbation matrix in harmonic oscillator basis
cmatrix_t *V_pert = cmatrix_alloc(N, N);
for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
        // Compute <i|x^4|j> using numerical integration
        double integral = 0.0;
        for (int k = 0; k < grid_N; k++) {
            double x4 = x[k]*x[k]*x[k]*x[k];
            integral += conj(psi_i[k]) * x4 * psi_j[k] * dx;
        }
        CMAT(V_pert, i, j) = c_real(lambda * integral);
    }
}

// Apply perturbation theory
perturb_result_t result = perturb_nondeg(eig0, n, V_pert, 1e-10);
```

### Degenerate Perturbation Theory

When two or more unperturbed states have the same energy, non-degenerate perturbation theory fails. The perturbation must be diagonalized within the degenerate subspace.

```c
eigen_t *perturb_degenerate(const double *energies, const cmatrix_t *V_pert,
                            const int *deg_indices, int deg_size) {
    // Build effective Hamiltonian in degenerate subspace
    cmatrix_t *H_eff = cmatrix_alloc(deg_size, deg_size);
    for (int i = 0; i < deg_size; i++) {
        for (int j = 0; j < deg_size; j++) {
            CMAT(H_eff, i, j) = CMAT(V_pert, deg_indices[i], deg_indices[j]);
        }
    }
    return cmatrix_eigh_generic(H_eff);
}
```

### Fermi's Golden Rule

For time-dependent perturbations, the transition rate from state _i_ to _f_ is:

$$
W_{i \rightarrow f}=\frac{2\pi}{\hbar}|\langle f∣V∣i \rangle∣^2 \rho(E_f)
$$

```c
double fermi_golden_rate(const cmatrix_t *V_pert, int i, int f, double rho_E) {
    complex_t Vfi = CMAT(V_pert, f, i);
    return 2.0 * M_PI * c_abs2(Vfi) * rho_E;
}
```
