/*
 * Two-site Hubbard model (Hubbard dimer), solved by exact diagonalisation
 *
 *   H = -t sum_{s} (c+_{1s} c_{2s} + h.c.) + U (n_{1up} n_{1dn} + n_{2up}
 * n_{2dn})
 *
 * 16-dimensional Fock space is built with Jordan-Wigner signs. Mode bits: 0 =
 * site 1 up, 1 = site 1 down, 2 = site 2 up, 3 = site 2 down; a basis state is
 * 4-bit mask of occupied modes
 *
 * At half filling (N = 2) spectrum is E = U/2 -+ sqrt(U^2/4 + 4 t^2) (two
 * singlets), U (antisymmetric charge-transfer singlet, which does not hop) and
 * 0 (triplet, three-fold), so ground state is a singlet that crosses over from
 * a band-like state (U = 0) to a Heisenberg pair with J = 4 t^2 / U (U >> t)
 * A state of definite N and S_z has a diagonal one-site reduced density matrix,
 * whose entropy measures charge and spin entanglement between two sites
 *
 * States are complex amplitude arrays re[16], im[16] indexed by mask
 */
#ifndef HUBBARD_DIMER_H
#define HUBBARD_DIMER_H

enum { HD_DIM = 16 };

typedef struct {
  double t, u;
  double e[HD_DIM];         /* energies, ascending */
  double v[HD_DIM][HD_DIM]; /* v[k][mask]: amplitude of basis state `mask` in
                               eigenstate k */
  int n[HD_DIM];            /* particle number of eigenstate k */
  int twosz[HD_DIM];        /* 2 S_z of eigenstate k */
} hd_solution_t;

typedef struct {
  double energy, double_occ, ss,
      entropy; /* <H>, <n1up n1dn + n2up n2dn>, <S1.S2>, thermodynamic entropy S
                  / k_B */
} hd_thermal_t;

/* Diagonalises H; 0 on success, -1 for non-finite parameters */
int hd_solve(double t, double u, hd_solution_t *s);

/* Fills idx with six N = 2 eigenstates in ascending energy order (idx[0] is
 * ground state) */
void hd_n2_states(const hd_solution_t *s, int idx[6]);

/* Basis state |mask> as (re, im). */
void hd_basis_state(int mask, double *re, double *im);

/* Eigenstate k as (re, im). */
void hd_eigenstate(const hd_solution_t *s, int k, double *re, double *im);

/* exp(-i H time) |psi0> (hbar = 1); -1 on non-finite time. */
int hd_evolve(const hd_solution_t *s, const double *re0, const double *im0,
              double time, double *re, double *im);

double hd_expect_energy(const hd_solution_t *s, const double *re,
                        const double *im);
double hd_expect_double_occ(const double *re,
                            const double *im); /* <n1up n1dn + n2up n2dn> */
double hd_expect_sz(const double *re, const double *im,
                    int site); /* <(n_up - n_dn) / 2> of site 0 or 1 */
double hd_expect_ss(const double *re, const double *im); /* <S1 . S2> */

/* Occupation probabilities of one site: p[0] empty, p[1] up, p[2] down, p[3]
 * double */
void hd_site_probs(const double *re, const double *im, int site, double p[4]);

/* Entropy -sum p ln p of those probabilities: entanglement entropy of a site
 * for a state of definite N and S_z */
double hd_site_entropy(const double *re, const double *im, int site);

/* Canonical (N = 2) thermal averages at temperature temp (k_B = 1); -1 unless
 * temp > 0 */
int hd_thermal_n2(const hd_solution_t *s, double temp, hd_thermal_t *out);

#endif
