#ifndef QMC_TFIM_QUENCH_H
#define QMC_TFIM_QUENCH_H

/*
 * Quantum quench of transverse-field Ising chain, solved exactly with free
 * fermions
 *
 *   H = -J \sum_j sz_j sz_(j+1) - h \sum_j sx_j        (ring of N sites)
 *
 * Chain starts in ground state (even parity) at field h_i and field is switched
 * suddenly to h_f. After a Jordan-Wigner transformation each pair of momenta
 * (k, -k), k = (2m + 1) pi / N (antiperiodic, m = 0 .. N/2 - 1), is a spin-1/2
 * with Bloch vector s_k. It starts anti-parallel to field direction n_k(h_i)
 * and precesses about n_k(h_f) at angular frequency 2 eps_k, with
 *
 *   n_k(h) = (0, J sin k, h - J cos k) / |...|
 *   \eps_k = 2 \sqrt(J^2 + h^2 - 2 J h \cos(k))
 *
 * (eps_k is quasiparticle energy of ising_exact_ground_energy_per_site).
 * Everything below is built from these Bloch vectors, so cost is
 * polynomial in N rather than 2^N:
 *
 *   transverse magnetisation   <sx>(t) = -(2/N) \sum_k s_z(k, t)
 *   Loschmidt echo             L(t) = prod_k [1 - (1 - (n_i . n_f)^2)
 *                              \sin^2(\eps_k t)]
 *   rate function              \lambda(t) = -ln L(t) / N
 *   correlations               <sz_0 sz_r>(t) = (-1)^r Pf(R)
 *
 * where R is 2r x 2r antisymmetric matrix of Majorana contractions built from
 * Fourier sums of s_k. When quench crosses critical point h = J rate function
 * has kinks at t*_n = \pi (n + 1/2) / \eps_k*, with cos k* = (J^2 + h_i h_f) /
 * (J (h_i + h_f)): dynamical quantum phase transitions. All results agree with
 * exact diagonalisation (ising_chain.h) at small N, which is how module is
 * tested
 */

enum { TFIM_QUENCH_MAX_SITES = 4096, TFIM_QUENCH_MAX_R = 64 };

/* Bloch vectors s_k(t) of N/2 momentum pairs, three numbers per mode
 * (x, y, z) in s[3 m .. 3 m + 2], for k = (2 m + 1) pi / N. Returns 0, or -1
 * for NULL s, odd or out-of-range n_sites (2 .. TFIM_QUENCH_MAX_SITES), or
 * non-finite J, h_i, h_f, t */
int tfim_quench_modes(double J, double h_i, double h_f, double t, int n_sites,
                      double *s);

/* Transverse magnetisation <sx> after quench (t = 0 gives ground state of h_i)
 * NaN on invalid input */
double tfim_quench_mx(double J, double h_i, double h_f, double t, int n_sites);

/* Loschmidt rate function -\ln L(t) / N, with L(t) =
 * |<\psi_i|\exp(-iHt)|\psi_i>|^2
 *NaN on invalid input */
double tfim_quench_loschmidt_rate(double J, double h_i, double h_f, double t,
                                  int n_sites);

/* Correlation <sz_0 sz_r>(t) for r = 0 .. r_max into zz (zz[0] = 1). r_max must
 * be at most min(n_sites / 2 - 1, TFIM_QUENCH_MAX_R) so string does not wrap
 * around ring. Returns 0, or -1 on invalid input, -2 on allocation failure */
int tfim_quench_zz(double J, double h_i, double h_f, double t, int n_sites,
                   int r_max, double *zz);

/* Critical time t*_n (n = 0, 1, ...) of dynamical quantum phase transition for
 * a quench across h = J in thermodynamic limit. NaN if quench does not cross
 * critical point, for n < 0 or non-finite input */
double tfim_quench_critical_time(double J, double h_i, double h_f, int n);

#endif
