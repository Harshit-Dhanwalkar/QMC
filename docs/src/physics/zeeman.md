# Zeeman Effect

The weak-field (anomalous) Zeeman effect for a single valence electron ($s=\frac12$ fixed), building on the [angular momentum coupling](angular_momentum.md) and [fine structure](fine_structure.md) machinery. Implemented in `physics/zeeman.h`, using the same doubled-quantum-number convention (`j_2 = 2j`, `mj_2 = 2*mj`).

## Physics

To first order in the external field $B$ (valid when the Zeeman splitting is small compared to the fine-structure splitting):

$$
H' = \mu_B B\,\frac{L_z + 2S_z}{\hbar} = \mu_B B\,\frac{J_z+S_z}{\hbar}
$$

Giving an energy shift $\delta E = g_J\,\mu_B B\,m_j$, where the Land$\'e$
g-factor is:

$$
g_J = 1 + \frac{j(j+1)+s(s+1)-l(l+1)}{2j(j+1)}, \qquad s=\tfrac12
$$

## Implementation

```c
/* Returns NAN if j_2 isn't an allowed coupling of l and spin-1/2
   (j_2 != 2l+1 and j_2 != 2l-1), or if j=0. */
double zeeman_lande_g_factor(int l, int j_2);

/* dE = g_J * mu_B * B * mj, mj = mj_2/2. mu_B caller-supplied
   (use 1.0 for natural units). */
double zeeman_energy_shift(int l, int j_2, int mj_2, double B, double mu_B);
```

```c
int j_2 = 2*1 + 1; // j = l + 1/2 branch, l=1 -> j=3/2 -> j_2=3
double g = zeeman_lande_g_factor(1, j_2);
double dE = zeeman_energy_shift(1, j_2, /*mj_2=*/1, B, 1.0);
```

### Independent Cross-Check

```c
/* <Sz>/hbar = sum_{ml,ms} |<l ml; s ms | j mj>|^2 * ms, computed directly
   from couple_states()'s Clebsch-Gordan expansion and not from g_J.
   Should equal (g_J - 1) * mj to numerical precision.
   Returns NAN if (l, j_2, mj_2) is an invalid coupling. */
double zeeman_sz_expect_from_coupling(int l, int j_2, int mj_2);
```

This mirrors the validation pattern already used in [Fine Structure](fine_structure.md)'s `spin_orbit_ls_expect_from_coupling` — an independent recomputation from the Clebsch-Gordan coupling coefficients, cross-checked against the closed-form result rather than trusting a single derivation. Since $\langle S_z\rangle/\hbar = (g_J-1)m_j$ follows algebraically from the same angular momentum algebra that gives $g_J$ itself, agreement between `zeeman_sz_expect_from_coupling` and `(zeeman_lande_g_factor(l, j_2) - 1) * mj` confirms both the `couple_states`/Clebsch-Gordan implementation and the closed-form Land$\'e$ formula are consistent with each other.

## Running the Example

```sh
./build/eg_26_zeeman
```
