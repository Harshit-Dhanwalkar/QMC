# Tight-Binding Models and Band Topology

Bloch Hamiltonians, band structures, density of states and topological invariants for multi-band models on 1D and 2D lattices. Natural units ($\hbar = 1$). Implemented in `physics/tight_binding.h`.

## The model

A unit cell holds `n_orb` orbitals, all placed at the cell origin. Each term added with `tb_add_hopping` is $t\,c^\dagger_{i,R}\,c_{j,R+n_1a_1+n_2a_2} + \text{h.c.}$, and each orbital has an on-site energy. The Bloch Hamiltonian

$$
H_{ij}(\mathbf k) \mathrel{+}= t\,e^{i(k_1n_1 + k_2n_2)}, \qquad H_{ji}(\mathbf k) \mathrel{+}= \overline{t\,e^{i(k_1n_1 + k_2n_2)}}
$$

is exactly periodic, $H(\mathbf k + 2\pi) = H(\mathbf k)$. Wave vectors are in reduced coordinates $(k_1, k_2) = (\mathbf k\cdot a_1, \mathbf k\cdot a_2)$, so the Brillouin zone is $[0, 2\pi)^2$ for every lattice and no geometry is needed.

## Topological invariants

The invariants use the gauge-invariant link-variable method of Fukui, Hatsugai and Suzuki: eigenvectors on a uniform $n_k$ grid are compared only through determinants of overlap matrices, so no smooth gauge is needed and several bands can be treated as a group.

- **Berry flux** of a plaquette is $F = \arg(U_1U_2U_3U_4)$, the product of four link variables around it; `tb_berry_flux_grid` returns the whole grid, which shows where the Berry curvature lives.
- **Chern number** $C = \frac{1}{2\pi}\sum F$ (`tb_chern_number`). It is returned as a `double` so its distance from an integer shows whether the grid has converged.
- **Wilson-loop (Zak) phase** $\gamma = -\arg\prod_k\det\langle u_k|u_{k+1}\rangle$ along $k_1$ (`tb_wilson_phase`), and its winding as $k_2$ runs round the zone (`tb_wilson_winding`), which equals the Chern number as an independent check.

A band group must be gapped from the rest of the spectrum; otherwise the functions return `TB_ERR_DEGENERATE`.

## Ready-made models

| Constructor                           | Model                                                                                                                     |
| ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `tb_model_chain(t)`                   | 1D chain, $E = 2t\cos k$                                                                                                  |
| `tb_model_ssh(v, w)`                  | Su–Schrieffer–Heeger chain; Zak phase $\pi$ when $\lvert w\rvert > \lvert v\rvert$                                        |
| `tb_model_square(t)`                  | square lattice, $E = 2t(\cos k_1 + \cos k_2)$                                                                             |
| `tb_model_graphene(t)`                | honeycomb lattice with Dirac points at $K$ and $K'$                                                                       |
| `tb_model_haldane(t1, t2, phi, mass)` | Haldane model; lower-band Chern number $\text{sign}(t_2\sin\phi)$ when $\lvert M\rvert < 3\sqrt3\lvert t_2\sin\phi\rvert$ |
| `tb_model_kagome(t)`                  | kagome lattice with an exactly flat band at $E = -2t$                                                                     |
| `tb_model_hofstadter(p, q, t)`        | square lattice with flux $p/q$ per plaquette as a $q$-orbital magnetic cell                                               |

Analytic references are provided for the dispersions (`tb_chain_energy`, `tb_square_energy`, `tb_graphene_energy`), the Haldane valley gap and Chern number (`tb_haldane_valley_gap`, `tb_haldane_chern`), the SSH Zak phase (`tb_ssh_zak_phase`) and the kagome flat band (`tb_kagome_flat_band_energy`).

## Hofstadter and the quantum Hall effect

With flux $p/q$ per plaquette the single square-lattice band splits into $q$ sub-bands. When the Fermi level sits in a gap, the Hall conductance is the integer $\sigma_{xy} = C\,e^2/h$, where $C$ is the sum of the Chern numbers of the bands below the gap (the TKNN result). Plotting every allowed energy against the flux gives the self-similar butterfly.

## Example

```c
tb_model_t *m = tb_model_hofstadter(1, 3, 1.0);   /* flux 1/3, three bands */
for (int band = 0; band < 3; band++) {
  double c;
  tb_chern_number(m, band, 1, 24, &c);
  printf("band %d: Chern = %+.6f\n", band, c);
}
tb_model_free(m);
```

The three bands carry Chern numbers $+1$, $-2$ and $+1$, which sum to zero as they must. For the SSH chain, `tb_wilson_phase` on the lower band gives $0$ for $v = 1, w = 0.5$ and $\pi$ for $v = 1, w = 2$. The full program is `examples/eg_75_tight_binding.c`.

## Verification

`tests/test_tight_binding.c` checks the chain, square and graphene dispersions against closed forms (with the graphene Dirac points, band lines and DOS), the Haldane valley gaps and Chern phase diagram, the SSH Zak phase against the edge states of an open chain, the kagome flat band, the Hofstadter Bloch matrix and band Chern numbers against `lattice.h`, and that the Wilson-loop winding equals the Chern number.
