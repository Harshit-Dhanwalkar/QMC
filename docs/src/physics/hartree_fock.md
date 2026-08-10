# Hartree-Fock (Restricted, s-Orbitals, Self-Consistent Field (SCF))

The Hartree-Fock method approximates the wave function of an $N$-body interacting fermion system using a single Slater determinant of spin-orbitals.

Restricted Hartree-Fock (RHF) self-consistent field for closed-shell atoms/ions whose occupied subshells are all s-type ($l=0$): helium ($1s^2$), beryllium ($1s^2 2s^2$), and similar s-shell ions. Atomic units throughout ($\hbar = m_e = e = 4\pi\varepsilon_0 = 1$). Implemented in `physics/hartree_fock.h`.

## Why s-orbitals only

For a central potential, the electron-electron interaction $1/|\mathbf{r}-\mathbf{r}'|$ expands in Legendre multipoles $Y^L$. Angular (Gaunt-coefficient) selection rules restrict which multipole $L$ actually contributes to the direct ($J$) and exchange ($K$) integrals between two orbitals of angular momenta $l_a, l_b$. For $l_a = l_b = 0$, the only surviving multipole is $L=0$, so both $J$ and $K$ collapse onto the _same_ radial kernel $1/\max(r,r')$ - this is what makes the exchange operator an exact, directly discretizable radial expression rather than requiring a general multipole expansion. That simplification is specific to s-orbitals; it's why this solver doesn't generalize to $p$, $d$, ... subshells as-is.

## Model

$n_{\text{orbitals}}$ distinct radial s-orbitals $u_1(r) < u_2(r) < \ldots$ (with $u = rR$, normalized $\int u_k(r)^2\,dr = 1$), each doubly occupied (spin up + spin down, Aufbau order) - i.e. this models atoms/ions with electron count $= 2 \times n_{\text{orbitals}}$, all in s subshells (He: `n_orbitals=1`, Be: `n_orbitals=2`, ...).

The closed-shell RHF Fock operator (Szabo & Ostlund eq. 3.184):

$$
F = h + \sum_k (2J_k - K_k), \qquad h = -\frac12\frac{d^2}{dr^2} + \frac{l(l+1)}{2r^2} - \frac{Z}{r}
$$

discretized as a dense real-symmetric matrix and diagonalized each SCF iteration via a Hermitian eigensolver. Total electronic energy:

$$
E = \sum_k 2\varepsilon_k - \sum_{i,j}(2J_{ij} - K_{ij})
$$

Where $\varepsilon_k$ are the converged orbital (Fock) eigenvalues.

## Implementation

```c
typedef struct {
  int n_orbitals;           // number of doubly-occupied s-orbitals
  int N;                    // radial grid size
  double Z;                 // nuclear charge used
  double *orbital_energies; // size n_orbitals, converged Fock eigenvalues
  cvector_t **orbitals;     // size n_orbitals; u_k(r) = r*R_k(r) in .re,
                            //  normalized: integral u_k(r)^2 dr = 1
  double total_energy;      // converged total electronic energy, Hartree
  int iterations;           // SCF iterations actually performed
  int converged;            // 1 if converged within max_iter, else 0
} hf_result_t;

hf_result_t *hartree_fock_atom_s_orbitals(double *r, int N, double Z,
                                          int n_orbitals, double mix,
                                          double tol, int max_iter);

void hf_result_free(hf_result_t *res);
```

`r` should be a **uniform** radial grid with `r[0] > 0` (avoiding the Coulomb singularity - same convention as [Central Potentials](central_potential.md)'s `central_potential_radial_solve`), and `N >= 10` or so for anything meaningful. `mix` is a linear density-mixing fraction in $(0,1]$ (`new = (1-mix)*old + mix*new_raw`) applied each SCF update for stability - `1.0` means no damping; smaller values (0.3–0.5) trade convergence speed for stability against oscillation.

```c
double *r = linspace(1e-4, 20.0, 2000);
hf_result_t *he = hartree_fock_atom_s_orbitals(r, 2000,
                                               /*Z=*/2.0,
                                               /*n_orbitals=*/1,
                                               /*mix=*/0.5,
                                               1e-8,
                                               200);
if (he->converged) {
    printf("Helium RHF energy: %.6f Hartree\n", he->total_energy);
}
hf_result_free(he);
```

Always check `hf_result_t.converged` before trusting `total_energy` - `max_iter` can be exhausted without convergence, in which case the returned struct still holds the last-iterate values.

## Running the Example

```sh
./build/eg_20_hartree_fock
```

## Relation to the simpler variational approach

[Helium & Two-Electron Atoms](helium.md) covers a much simpler effective-nuclear-charge variational treatment restricted to two electrons with a fixed single-parameter trial wavefunction. This module generalizes that in two ways at once: an arbitrary number of doubly-occupied s subshells (not just $1s^2$), and self-consistent orbital shapes solved numerically rather than a single closed-form charge parameter - at the cost of losing the closed-form solution [Helium & Two-Electron Atoms](helium.md)'s variational method has.
