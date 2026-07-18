#ifndef QMC_HELIUM_H
#define QMC_HELIUM_H

/*
 * Ground-state energy of a two-electron atom/ion via classic
 * effective-nuclear-charge variational method (Griffiths "Introduction to
 * Quantum Mechanics"; Bransden & Joachain Ch. 5-7).
 *
 * HACK: Trial wavefunction: \psi(r1,r2) = \phi(r1) * \phi(r2), \phi(r;Z')
 * hydrogen-like 1s orbital with effective charge Z' (spatially  symmetric,
 * appropriate for the spin-singlet ground state - electron spin/antisymmetry is
 * not modeled explicitly here, only spatial variational energy).
 *
 * All energies in Hartree atomic units (\hbar = m_e = e = 4* pi * EPS0 = 1).
 * convert Hartree -> eV
 *
 * Using known hydrogenic-1s expectation values for TRIAL charge Z'
 * (via the virial theorem for effective Hamiltonian
 * H'=-1/2 grad^2 - Z'/r that \phi(r;Z') is an eigenfunction of):
 *   <T> = Z'^2/2 per electron   (kinetic operator is charge-independent)
 *   <1/r> = Z'                  (so <-Z/r> = -Z*Z' with real charge Z)
 *   <1/r12> = (5/8)*Z'          (standard two-center Coulomb integral)
 * gives closed-form total energy:
 *   E(Z') = Z'^2 - 2*Z*Z' + (5/8)*Z'
 * which is exactly minimized at Z'_opt = Z - 5/16, giving
 *   E_opt = -(Z - 5/16)^2
 * For helium (Z=2): Z'_opt=1.6875, E_opt=-2.84765625 Hartree
 * (~-77.5 eV),
 * HACK: variational theorem guarantees E_opt >= E_true, and gap here is known
 * cost of neglecting electron correlation beyond this simple product-orbital
 * ansatz.
 */

// E(Z') in Hartree, for nuclear charge Z (2 for neutral helium).
double helium_variational_energy(double Z_eff, double Z);

// Exact analytic minimizer: Z'_opt = Z - 5/16.
double helium_optimal_zeff_analytic(double Z);

// Exact analytic minimum energy: E_opt = -(Z-5/16)^2, in Hartree.
double helium_ground_state_energy_analytic(double Z);

/*
 * Numerically minimizes helium_variational_energy over Z' using
 * generic golden_section_minimize
 * Search bounds are [max(0.01, 0.1*Z), 2*Z], which comfortably brackets
 * Z'_opt=Z-5/16 for any physically sensible Z >= 1.
 *
 * If Zeff_opt_out is non-NULL, optimal Z' found is written there.
 * Returns the minimized energy in Hartree.
 */
double helium_ground_state_energy_numeric(double Z, double tol,
                                          double *Zeff_opt_out);

#endif
