# Density Matrix Renormalization Group (DMRG)

DMRG variationally approximates the ground state of a 1D quantum chain by
growing and repeatedly truncating a real-space block basis, keeping only the
$m$ most-probable states of the block's reduced density matrix at each step.
Implemented for the open-boundary spin-1/2 XXZ chain:

$$
H = \sum_{j=0}^{N-2} \Big[ J_z S^z_j S^z_{j+1} + \frac{J_{xy}}{2}\big(S^+_j S^-_{j+1} + S^-_j S^+_{j+1}\big) \Big]
$$

in `physics/dmrg.h` (infinite-system algorithm, White, Phys. Rev. Lett. 69,
2863 (1992)) and `physics/finite_dmrg.h` (finite-size sweeps on top of it,
White, Phys. Rev. B 48, 10345 (1993)).

## Infinite-system algorithm

Starting from a single-site block (dim 2, $H=0$):

1. **Enlarge**: new block = old block $\otimes$ one new site, with the
   coupling bond between them folded into the new block Hamiltonian.
2. **Superblock**: form enlarged block $\otimes$ its own mirror image
   (reflection symmetry stands in for an explicit environment block, since
   there's no separate, already-built environment yet) and find its ground
   state by dense diagonalization.
3. **Truncate**: diagonalize the ground state's reduced density matrix on
   the enlarged-block subspace; keep the $m$ eigenvectors of largest
   eigenvalue. This is the density-matrix truncation the method is named
   for.
4. Repeat from step 1 until the superblock's total site count reaches the
   target $N$ (each step grows the chain by 2 sites, one per mirrored side,
   so $N$ is reached exactly when even, or overshot by 1 when odd).

For a gapped chain this reaches machine precision once $m$ is large enough.
For a fixed, insufficient $m$ on a critical/gapless chain (e.g. isotropic
Heisenberg, $J_z=J_{xy}$), every intermediate block was built against an
environment that was itself still small and approximate this is exactly
what finite-size sweeping (below) improves on.

```c
dmrg_result_t *dmrg_run(int N_target, double Jz, double Jxy, int m_max);
```

`N_target` is rounded up internally to the nearest even number. Returns
`NULL` for invalid input (`N_target < 2`, `m_max < 1`, non-finite `Jz`/`Jxy`)
or allocation failure; free a successful result with plain `free()`.

```c
dmrg_result_t *r = dmrg_run(/*N_target=*/40, /*Jz=*/1.0, /*Jxy=*/1.0,
                            /*m_max=*/20);
printf("E/N = %.8f (truncation error %.2e)\n", r->energy_per_site,
      r->truncation_error);
free(r);
```

## Finite-size sweeps

The infinite algorithm never revisits a block once built. Finite-size
sweeping fixes the total chain length $N$ up front, keeps a cache of every
block size $1 \ldots N-1$ seen so far (reflection symmetry means one cache
serves as both left- and right-block storage), and repeatedly
re-diagonalizes superblocks built from one *growing* side plus a
previously-cached *environment* side of the complementary size,
re-truncating only the growing side each step:

1. **Warmup**: run the infinite algorithm to build cache entries for every
   size $1 \ldots N/2$, truncating each one (including the last, which
   `dmrg_run` itself leaves untouched since it never needs it again - here
   it must be capped like every other entry or the first sweep's
   convergence checkpoint looks artificially better than it is).
2. **Sweep**: pick a growing size $l$ starting at 2, form the superblock
   from `enlarge(cache[l-1])` (system) and `cache[N-l]` (environment, read
   only), find its ground state, truncate the enlarged side to at most
   `m_max` states, store the result back into `cache[l]`. Advance $l$ to
   $N-1$.
3. **Reverse**: the shrunk side now becomes the side that grows, using the
   fresh blocks the other side just produced as its improving environment;
   sweep back down to $l=1$. One outward-and-back pass is one full sweep.
4. Report the ground energy at the most balanced bipartition ($l=N/2$)
   after each full sweep.

Because every sweep step is itself an exact diagonalization of the current
superblock, the reported central-bond energy is non-increasing (up to
solver round-off) as more sweeps are performed, and matches or improves on
the infinite algorithm's result for the same $(N, m)$ once several sweeps
have run.

```c
typedef struct {
  double energy;            // ground energy at final l=N/2 bipartition
  double energy_per_site;
  int N_reached;            // total chain length actually used (even)
  int m_max;
  int n_sweeps;              // full sweeps actually performed
  double truncation_error;   // discarded RDM weight on the final l=N/2 step
  double *sweep_energy;      // length n_sweeps: l=N/2 energy after each sweep
} finite_dmrg_result_t;

finite_dmrg_result_t *finite_dmrg_run(int N_target, double Jz, double Jxy,
                                      int m_max, int n_sweeps);
void finite_dmrg_result_free(finite_dmrg_result_t *r);
```

```c
finite_dmrg_result_t *r = finite_dmrg_run(40, 1.0, 1.0, /*m_max=*/8,
                                          /*n_sweeps=*/5);
for (int s = 0; s < r->n_sweeps; s++) {
    printf("sweep %d: E = %.10f\n", s + 1, r->sweep_energy[s]);
}
finite_dmrg_result_free(r);
```

`N_target` is rounded up to the nearest even number, same convention as
`dmrg_run`. Returns `NULL` for invalid input (`N_target < 4`, `m_max < 1`,
`n_sweeps < 1`, non-finite `Jz`/`Jxy`) or allocation failure.

## Running the Examples

```sh
./build/eg_69_dmrg
./build/eg_70_finite_dmrg
```

The second example plots central-bond energy against sweep number (should
decrease monotonically then plateau) and compares finite-sweep vs.
infinite-algorithm accuracy at fixed, small `m` across a range of chain
lengths.

