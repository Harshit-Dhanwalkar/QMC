# Coupled Cluster Singles and Doubles (CCSD)

Post-Hartree-Fock correlation method that exponentiates single and double
excitation operators, $|\Psi_{\text{CCSD}}\rangle = e^{T_1+T_2}|\Phi_0\rangle$,
acting on a converged closed-shell RHF reference. More accurate (and more
expensive: $O(n^6)$ per iteration) than MP2 or CISD, and size-extensive unlike
truncated CI. Implemented in `physics/ccsd.h`/`physics/ccsd.c`, on top of
`molecular_hf.c`'s RHF and `molecular_ao_to_mo`'s MO-basis integrals.

## Formulation

Uses the standard spin-orbital formulation (Stanton, Gauss, Watts & Bartlett,
J. Chem. Phys. 94, 4334 (1991)): antisymmetrized spin-orbital integrals are
built from the spatial RHF MOs (doubled into alternating alpha/beta spin
orbitals, $n_{so} = 2 n_{\text{spatial}}$), and the $T_1$/$T_2$ amplitude
equations are solved by fixed-point iteration through the usual
$F_{ae}, F_{mi}, F_{me}, W_{mnij}, W_{abef}, W_{mbej}$ intermediates.

This is not spin-adapted (working directly in the smaller spatial-orbital
basis, as production codes do): the equations are simpler to state and
verify correctly at the cost of $n_{so}^4 = 16\, n_{\text{spatial}}^4$
storage/compute scaling instead of $n_{\text{spatial}}^4$. A full
spin-adapted rewrite remains a possible, much larger, future project.

Within the spin-orbital formulation, three of the original seven internal
$n_{so}^4$ arrays are eliminated without changing the arithmetic at all:
the antisymmetrized integral tensor $V$ and the trivial energy-denominator
tensor $D_{ijab}$ are computed on demand instead of cached (`v_elem()` /
`d_ijab()` in `ccsd.c`), and $W_{abef}$, the single largest intermediate
whenever $n_v > n_o$, is fused directly into its one point of use instead
of being built and stored every iteration. The remaining per-iteration cost
is spread across cores with OpenMP; because every remaining intermediate's
outer loop is embarrassingly parallel (each element computed start-to-finish
within one thread), results are bit-identical regardless of thread count.

## Implementation

```c
typedef struct {
  double correlation_energy; // E_CCSD - E_RHF
  double total_energy;       // E_RHF + correlation_energy
  int converged;
  int iterations;
} ccsd_result_t;

ccsd_result_t *ccsd_run(int n_spatial, const double *h_mo,
                        const double *eri_mo, const double *mo_energy,
                        int n_electrons, int n_frozen_spatial, double e_rhf,
                        double conv_tol, int max_iter);
```

`h_mo` is the MO-basis core Hamiltonian and `eri_mo` the $n_{\text{spatial}}^4$
chemist-notation ERI tensor, both from `molecular_ao_to_mo`; `mo_energy` and
`e_rhf` come from the converged `molecular_hf_result_t`. `n_frozen_spatial`
excludes that many lowest-energy spatial orbitals from correlation (0 for
none; must satisfy `2*n_frozen_spatial < n_electrons`). Returns `NULL` for
invalid input (odd electron count, frozen space too large, non-positive
`n_spatial`) or allocation failure.

```c
molecular_hf_result_t *rhf = molecular_hf_run(/* ... */);
double *h_mo, *eri_mo; // from molecular_ao_to_mo, using rhf->C
ccsd_result_t *cc = ccsd_run(n_spatial, h_mo, eri_mo, rhf->orbital_energies,
                             n_electrons, /*n_frozen_spatial=*/0,
                             rhf->total_energy, 1e-9, 100);
if (cc->converged) {
    printf("CCSD total energy: %.8f Hartree\n", cc->total_energy);
}
free(cc);
```

For downstream consumers that need the converged amplitudes themselves -
the perturbative $(T)$ triples correction (CCSD(T)) and
UCC - use `ccsd_run_ex` instead, which behaves identically to
`ccsd_run` when `amplitudes_out` is `NULL`:

```c
ccsd_amplitudes_t *amp = NULL;
ccsd_result_t *cc = ccsd_run_ex(n_spatial, h_mo, eri_mo,
                                rhf->orbital_energies, n_electrons, 0,
                                rhf->total_energy, 1e-9, 100, &amp);
// amp->V is only materialized here, once, because it was asked for -
// ccsd_run's plain energy-only path never builds it at all.
ccsd_amplitudes_free(amp);
free(cc);
```

## Running the Example

```sh
./build/eg_50_ccsd
```
