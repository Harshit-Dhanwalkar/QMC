#ifndef QMC_MOLECULAR_INTEGRALS_H
#define QMC_MOLECULAR_INTEGRALS_H

/*
 * General Gaussian-type-orbital (GTO) molecular integrals via
 * McMurchie-Davidson (MD) Hermite-expansion scheme (McMurchie & Davidson, J.
 * Comput. Phys. 26, 218 (1974); see also Helgaker, Jorgensen & Olsen,
 * "Molecular Electronic-Structure Theory", ch. 9).
 * NOTE: This is lets second_quant.c / vqe.c build a real molecular
 * electronic-structure Hamiltonian (e.g. for H2 VQE) instead of only the toy
 * fermionic-lattice Hamiltonians second_quant_build_hopping_hamiltonian()
 * provides.
 *
 * Supports arbitrary Cartesian angular momentum (s, p, d, f, ...) and arbitrary
 * contraction (number of primitives), so this is a general basis-set /
 * general-molecule engine, not restricted to any one element or minimal basis.
 * WARN: though only basis-set builder shipped here (molint_basis_sto3g_h) is
 * STO-3G hydrogen, since that is exactly what's needed for a first real H2 VQE
 * target. Other basis functions can be constructed directly via
 * basis_function_t without any new engine.
 *
 * Units: atomic units throughout (bohr for length, Hartree for energy)
 */

#include "../core/matrix.h"

/* Boys function F_n(x) = integral_0^1 t^{2n} * \exp(-x t^2) dt, universal
 * building block of every GTO-based Coulomb-type integral (nuclear attraction,
 * ERIs). Computed via the \exp(-x)*convergent-power-series form (all terms same
 * sign, no cancellation) for F_nmax, then the stable downward recursion
 *   F_{n-1}(x) = (2x * F_n(x) + \exp(-x)) / (2n-1) for n < nmax
 * (WARN: upward recursion is numerically unstable and must never be used)
 *
 * Where
 *  boys_function:       single value F_n(x).
 *  boys_function_array: fills F[0..nmax] in one downward-recursion pass.
 */
double boys_function(int n, double x);
void boys_function_array(int nmax, double x, double *F /* length nmax+1 */);

/* A single primitive Cartesian Gaussian, unnormalized:
 *   g(r) = x^l y^m z^n \exp(-\alpha |r - center|^2)
 * l+m+n is this primitive's total angular momentum (0=s, 1=p, 2=d, ...).
 * Individual (l,m,n) triples distinguish e.g. px/py/pz (1,0,0)/(0,1,0)/
 * (0,0,1), or the six Cartesian d components.
 */
typedef struct {
  int l, m, n;
  double alpha;
} gto_primitive_t;

/* A contracted basis function: a fixed linear combination of primitives sharing
 * one center and one (l,m,n). coefficients[] are expected to already include
 * each primitive's own normalization and the overall contraction normalization.
 */
typedef struct {
  int l, m, n;
  double center[3];
  int n_primitives;
  double *exponents;    /* length n_primitives */
  double *coefficients; /* length n_primitives */
} basis_function_t;

/* Allocate a basis function and copy in exponents/raw (unnormalized)
 * literature contraction coefficients
 *
 * Returns NULL on allocation failure or n_primitives <= 0.
 */
basis_function_t *basis_function_alloc(int l, int m, int n,
                                       const double center[3], int n_primitives,
                                       const double *exponents,
                                       const double *raw_coefficients);
void basis_function_free(basis_function_t *bf);

/* In place:
 *   - folds each primitive's own Cartesian-Gaussian normalization constant into
 *     coefficients[i],
 *   - rescales the whole contraction so self-overlap integral(bf * bf) = 1.
 */
void molint_normalize_contraction(basis_function_t *bf);

/* A molecule: point nuclear charges (atomic units, so charge is in units of
 * |e|, e.g. 1.0 for hydrogen) at fixed 3D centers. Same molecule_t can be
 * paired with any basis set.
 */
typedef struct {
  int n_atoms;
  double *charge;      /* length n_atoms */
  double (*center)[3]; /* length n_atoms */
} molecule_t;

molecule_t *molecule_alloc(int n_atoms, const double *charge,
                           const double center[][3]);
void molecule_free(molecule_t *mol);
/* Classical nuclear-nuclear repulsion energy :
 *  \Sum_{A<B} Z_A Z_B / |R_A - R_B|
 * piece of total electronic-structure energy that has nothing to do with basis
 * set or electrons.
 */
double molecule_nuclear_repulsion(const molecule_t *mol);

/*
 * Core two-center / three-center / four-center integrals over a pair (or
 * quadruple) of contracted basis_function_t, using McMurchie-Davidson expansion
 * over every primitive combination in contraction(s). All take fully general
 * (l,m,n) on every center : s/p/d/... all go through the same code path, there
 * is no separate "s-function fast path".
 */

// <a|b>
double gto_overlap(const basis_function_t *a, const basis_function_t *b);

// <a|-1/2 nabla^2|b>
double gto_kinetic(const basis_function_t *a, const basis_function_t *b);

// <a|1/|r-C||b>
/* NOTE: Positive convention :integral itself, not yet multiplied by a nuclear
 * charge or a minus sign; molecular_nuclear_matrix below applies -Z_A per atom
 * and sums over atoms, matching one-electron Hamiltonian h = T - \sum_A Z_A *
 * this_integral).
 */
double gto_nuclear_attraction(const basis_function_t *a,
                              const basis_function_t *b,
                              const double center[3]);

/* (ab|cd) in chemists' notation: integral integral a(1)b(1) (1/r12)
 * c(2)d(2) d1 d2.
 * Full four-center, four-arbitrary-angular-momentum two-electron repulsion
 * integral.
 */
double gto_eri(const basis_function_t *a, const basis_function_t *b,
               const basis_function_t *c, const basis_function_t *d);

/*
 * Whole-basis matrix/tensor builders: loop the above over every pair (or
 * quadruple) of functions in a basis_function_t** array of length n_basis.
 *
 * Returned cmatrix_t are real-valued (stored as complex_t with zero imaginary
 * part).
 */
cmatrix_t *molecular_overlap_matrix(basis_function_t **basis, int n_basis);
cmatrix_t *molecular_kinetic_matrix(basis_function_t **basis, int n_basis);
/* h_core = T - \sum_A Z_A * gto_nuclear_attraction(., ., R_A),
 * NOTE: full one-electron Hamiltonian matrix (kinetic + nuclear attraction to
 * every atom in mol, not just one center). */
cmatrix_t *molecular_core_hamiltonian(basis_function_t **basis, int n_basis,
                                      const molecule_t *mol);

/* Flat row-major n_basis^4 array of (ij|kl) chemists'-notation ERIs,
 * index (((i * n_basis + j) * n_basis + k) * n_basis + l).
 * NOTE: Exploits 8-fold permutational symmetry (ij|kl) = (ji|kl) = (ij|lk) =
 * (kl|ij) = ... by computing each symmetry-distinct integral once and copying
 * to every equivalent slot, rather than calling gto_eri n_basis^4 times
 * independently.
 *
 * Returns NULL on n_basis<=0 or allocation failure. */
double *molecular_eri_tensor(basis_function_t **basis, int n_basis);
#define MOLINT_ERI(eri, n, i, j, k, l)                                         \
  (eri)[(((size_t)(i) * (n) + (j)) * (n) + (k)) * (n) + (l)]

/*
 * Convenience builder: STO-3G hydrogen 1s (Szabo & Ostlund Table 3.7 /
 * published EMSL/Gaussian94 STO-3G parameters), already normalized.
 *
 * Returns NULL on allocation failure. */
basis_function_t *molint_basis_sto3g_h(const double center[3]);

#endif
