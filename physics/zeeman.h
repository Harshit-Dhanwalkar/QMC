#ifndef QMC_ZEEMAN_H
#define QMC_ZEEMAN_H

/*
 * Weak-field (anomalous) Zeeman effect for a single valence electron
 * (s=1/2 fixed), in doubled-quantum-number convention j_2 = 2j, mj_2 = 2*mj.
 *
 * H' = \mu_B * B * (L_z + 2 S_z) / \hbar = \mu_B * B * (J_z + S_z) / \hbar,
 * to first order in B
 * (NOTE: valid when the Zeeman splitting is small compared to fine-structure
 * splitting).
 *
 * Energy shift: dE = g_J * \mu_B * B * mj
 * Lande g-factor: g_J = 1 + [j(j+1)+s(s+1)-l(l+1)] / (2 j(j+1)), s=1/2.
 */

/*
 * Lande g-factor for orbital l, total angular momentum j_2 (doubled),
 * s=1/2 fixed. Returns NAN if j_2 isn't an allowed coupling of l and
 * spin-1/2 (i.e. j_2 != 2l+1 and j_2 != 2l-1), or if j=0.
 */
double zeeman_lande_g_factor(int l, int j_2);

/*
 * Weak-field Zeeman energy shift: dE = g_J * mu_B * B * mj, mj = mj_2/2.
 * mu_B is caller-supplied; NOTE: 1.0 for natural units.
 */
double zeeman_energy_shift(int l, int j_2, int mj_2, double B, double mu_B);

/*
 * Independent cross-check of <S_z>/hbar for the coupled state
 * |l,s=1/2;j,mj>, computed directly from couple_states()'s Clebsch-
 * Gordan expansion into the uncoupled product basis |l,ml>|s,ms>:
 *   <Sz> = sum_{ml,ms} |<l ml; s ms | j mj>|^2 * ms
 * rather than from the closed-form g_J. Should equal (g_J - 1) * mj to
 * numerical precision
 *
 * Returns NAN if (l, j_2, mj_2) is an invalid coupling.
 */
double zeeman_sz_expect_from_coupling(int l, int j_2, int mj_2);

#endif
