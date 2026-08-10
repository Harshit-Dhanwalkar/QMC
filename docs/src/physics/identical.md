# Identical Particles

System of $N$ identical particles requires exchange symmetry under particle permutation operator $P_{ij}$.

Fermions and bosons require (anti)symmetrized many-body wavefunctions. `physics/identical.h` builds the machinery to evaluate an $N$-particle Slater determinant (fermions) or permanent (bosons) at a specific particle configuration.

## Exchange Symmetry

- **Bosons** (Symmetric states): $\psi(1, 2) = +\psi(2, 1)$
- **Fermions** (Antisymmetric states): $\psi(1, 2) = -\psi(2, 1)$

For two non-interacting particles in spatial single-particle states $\psi_a$ and $\psi_b$:

$$
\psi_{\pm}(\mathbf{r}_1, \mathbf{r}_2) = \frac{1}{\sqrt{2(1 \pm |\langle a|b \rangle|^2)}} \Big[ \psi_a(\mathbf{r}_1)\psi_b(\mathbf{r}_2) \pm \psi_b(\mathbf{r}_1)\psi_a(\mathbf{r}_2) \Big]
$$

## Exchange Energy

For spatial interactions $V(\mathbf{r}_1, \mathbf{r}_2)$:

$$
\langle V \rangle_{\pm} = J \pm K
$$

where:

- **Direct Integral**: $J = \int |\psi_a(\mathbf{r}_1)|^2 V(\mathbf{r}_1, \mathbf{r}_2) |\psi_b(\mathbf{r}_2)|^2 d^3r_1 d^3r_2$
- **Exchange Integral**: $K = \int \psi_a^*(\mathbf{r}_1) \psi_b^*(\mathbf{r}_2) V(\mathbf{r}_1, \mathbf{r}_2) \psi_b(\mathbf{r}_1) \psi_a(\mathbf{r}_2) d^3r_1 d^3r_2$

## Setup

Both constructions start from the same $N \times N$ matrix of single-particle orbitals sampled at each particle's position:

```c
// M_ij = orbitals[i]->data[indices[j]], size N x N
cmatrix_t *slater_matrix(cvector_t **orbitals, int N, const int *indices);
```

`orbitals[i]` is single-particle orbital $i$, tabulated on a common position grid of length $M$. `indices[j]` is the grid index of particle $j$'s position ($0 \le \texttt{indices[j]} < M$). Two particles at the same grid index, or two identical orbitals, both make the Slater determinant vanish exactly - this _is_ the Pauli exclusion principle falling out of the determinant structure, not a separately-enforced rule.

## Fermions: Slater determinant

$$
\psi(x_1,\ldots,x_N) = \frac{1}{\sqrt{N!}}\det\left[\phi_i(x_j)\right]
$$

```c
complex_t slater_determinant_value(cvector_t **orbitals, int N,
                                   const int *indices);
```

Implemented via $O(N^3)$ Gaussian elimination on the matrix from `slater_matrix`. Antisymmetric under exchange of any two particles by construction (determinant sign flips under row swap).

## Bosons: permanent

$$
\psi(x_1,\ldots,x_N) = \frac{1}{\sqrt{N!\prod_k n_k!}}\,\text{perm}\left[\phi_i(x_j)\right]
$$

```c
complex_t bosonic_permanent_value(cvector_t **orbitals, int N,
                                  const int *indices);
```

> Symmetric (invariant) under exchanging any two particles. $n_k$ is the occupation number of each distinct orbital among `orbitals[0..N-1]`; the $\prod_k n_k!$ factor only matters when multiple bosons share the same orbital - with all-distinct orbitals this reduces to the simpler $1/\sqrt{N!}$ normalization. Implemented via Ryser's formula, $O(2^N N)$ rather than the $O(N!)$ of a naive permanent expansion - still exponential, so this is only practical for modest $N$. `eg_25_boson_sampling.c` and `physics/boson_sampling.h` appear to build on this for a boson-sampling demo, but that header hasn't been reviewed yet - worth a page of its own in a later batch.

## Running the Example

```sh
./build/eg_09_identical_particles
```
