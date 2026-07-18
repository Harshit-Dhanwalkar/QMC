#ifndef QMC_FINE_STRUCTURE_H
#define QMC_FINE_STRUCTURE_H

/*
 * Hydrogen fine structure: spin-orbit coupling
 * closed-form total fine-structure formula for cross-checking.
 *
 * Conventions: l is the orbital quantum number (integer, l >= 0). Spin
 * is fixed at s=1/2 (electron). j_2 = 2*j and M_2 = 2*M, matching
 * angular.c's doubled-quantum-number convention (j=l+1/2 -> j_2=2l+1,
 * j=l-1/2 -> j_2=2l-1).
 */

// <L.S>/\hbar^2 = 0.5*(j(j+1) - l(l+1) - s(s+1)), s=1/2 fixed. Closed form.
double spin_orbit_ls_expect(int l, int j_2);

/*
 * <1/r^3>_{nl} for hydrogen = 1/(a0^3 n^3 l(l+1/2)(l+1)).
 * Returns 0 for l=0
 */
double hydrogen_expect_inv_r3(int n, int l, double hbar, double mass,
                              double e_charge, double eps0);

/*
 * Perturbative spin-orbit energy shift for hydrogen:
 *   dE_SO = \exp^2/(8 \pi \eps0 m^2 c^2) * <1/r^3>_{nl} * \hbar^2 *
 * <L.S>/\hbar^2 Returns 0 for l=0.
 *
 * NOTE: spin-orbit contribution ONLY. Total physical fine-structure shift also
 * includes relativistic kinetic-energy correction and (for l=0) Darwin term
 */
double hydrogen_spin_orbit_energy(int n, int l, int j_2, double hbar,
                                  double mass, double e_charge, double eps0,
                                  double c);

/*
 * Known closed-form TOTAL fine-structure shift:
 *   dE_fs = E_n * (alpha^2/n^2) * (n/(j+1/2) - 3/4)
 * where E_n = hydrogen_energy_level(n)
 */
double hydrogen_fine_structure_shift(int n, int j_2, double hbar, double mass,
                                     double e_charge, double eps0, double c);

/*
 * Direct verification of <L.S> using couple_states()'s CG coefficients,
 * evaluated via L.S = Lz Sz + 1/2*(L+ S- + L- S+) acting term-by-term on
 * uncoupled product basis. This does NOT call the closed-form
 * spin_orbit_ls_expect() - it recomputes <L.S> from scratch through
 * coupling coefficients and l_plus_op/l_minus_op, as independent
 * cross-check that couple_states() and ladder operators agree with
 * analytic j(j+1)-l(l+1)-s(s+1) result.
 *
 * Returns <L.S>/hbar^2 (dimensionless), or NAN if (l,j_2,M_2) is invalid.
 */
double spin_orbit_ls_expect_from_coupling(int l, int j_2, int M_2);

#endif
