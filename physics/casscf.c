/*
CASSCF (Complete Active Space SCF)
*/

#include "casscf.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "complex.h"
#include "second_quant.h"
#include <math.h>
#include <stdlib.h>

/* Direct (bit-manipulation) fermionic operators: To build RDMs directly from a
 * CI vector without constructing full operator matrices. MSB-first bit
 * convention: mode 0 is the leftmost/most significant bit of an n_modes-bit
 * Fock-space index
 */
static int casscf_count_bits_before(int state, int mode, int n_modes) {
  int count = 0;

  for (int k = 0; k < mode; k++) {
    int bitpos = n_modes - 1 - k;

    if ((state >> bitpos) & 1) {
      count++;
    }
  }

  return count;
}

static int casscf_direct_annihilate(int state, int mode, int n_modes,
                                    int *sign_out) {
  int bitpos = n_modes - 1 - mode;

  if (!((state >> bitpos) & 1)) {
    return -1;
  }

  int n_before = casscf_count_bits_before(state, mode, n_modes);
  *sign_out = (n_before % 2 == 0) ? 1 : -1;

  return state & ~(1 << bitpos);
}

static int casscf_direct_create(int state, int mode, int n_modes,
                                int *sign_out) {
  int bitpos = n_modes - 1 - mode;

  if ((state >> bitpos) & 1) {
    return -1;
  }

  int n_before = casscf_count_bits_before(state, mode, n_modes);
  *sign_out = (n_before % 2 == 0) ? 1 : -1;

  return state | (1 << bitpos);
}

/* Small dense linear algebra helpers (matrix exponential via
 * scaling-and-squaring + Taylor series). cmatrix_t is used throughout with zero
 * imaginary parts
 */
static cmatrix_t *casscf_identity(int n) {
  cmatrix_t *I = cmatrix_alloc(n, n);
  if (!I) {
    return NULL;
  }

  for (int i = 0; i < n * n; i++) {
    I->data[i] = c_zero();
  }
  for (int i = 0; i < n; i++) {
    CMAT(I, i, i) = c_real(1.0);
  }

  return I;
}

static double casscf_matrix_max_abs(const cmatrix_t *A) {
  double m = 0.0;

  for (int i = 0; i < A->nrows * A->ncols; i++) {
    double v = fabs(A->data[i].re);

    if (v > m) {
      m = v;
    }
  }

  return m;
}

/* \exp(K) for a small real matrix K, via scaling-and-squaring:
 *    \exp(K) = (\exp(K / 2^s)^(2^s)
 * with the inner exponential evaluated by a truncated Taylor series once
 * ||K/2^s|| is small enough for the series to converge rapidly. andard
 */
static cmatrix_t *casscf_matrix_exp(const cmatrix_t *K) {
  int n = K->nrows;
  double norm = casscf_matrix_max_abs(K);
  int s = 0;
  double scale = 1.0;

  while (norm * scale > 0.25 && s < 30) {
    scale *= 0.5;
    s++;
  }

  cmatrix_t *Ks = cmatrix_copy(K);
  cmatrix_scale(Ks, c_real(scale));

  cmatrix_t *result = casscf_identity(n);
  cmatrix_t *term = casscf_identity(n);

  for (int k = 1; k <= 15; k++) {
    cmatrix_t *new_term = cmatrix_multiply(term, Ks);

    cmatrix_free(term);
    term = new_term;
    cmatrix_scale(term, c_real(1.0 / k));

    cmatrix_t *new_result = cmatrix_add(result, term);

    cmatrix_free(result);
    result = new_result;
  }

  cmatrix_free(term);
  cmatrix_free(Ks);

  for (int i = 0; i < s; i++) {
    cmatrix_t *sq = cmatrix_multiply(result, result);

    cmatrix_free(result);
    result = sq;
  }

  return result;
}

/*
 * Active-space CI + RDM construction
 */
typedef struct {
  double energy; /* active-space CI ground energy, INCLUDING nuclear
                  * repulsion + frozen-core contribution */
  double *D1;    /* n_active x n_active */
  double *D2;    /* n_active^4, chemist-notation-compatible: D2[P,Q,R,S]
                  * pairs with eri_active[P,Q,R,S] in
                  * E = \sum h_eff * D1 + 1/2 sum eri * D2 */
} casscf_active_ci_t;

// Slice out the top-left n_sub x n_sub / n_sub^4 block of an n_basis-sized MO
// integral array (canonical ordering: core+active orbitals first, virtuals
// last, exactly as CISD/FCI's frozen-core mode already assumes)
static void casscf_slice_integrals(int n_basis, int n_sub, const double *h_mo,
                                   const double *eri_mo, double *h_sub,
                                   double *eri_sub) {
#define EMO_FULL(p, q, r, s)                                                   \
  eri_mo[(((size_t)(p) * n_basis + (q)) * (size_t)n_basis + (r)) *             \
             (size_t)n_basis +                                                 \
         (s)]

  for (int p = 0; p < n_sub; p++) {
    for (int q = 0; q < n_sub; q++) {
      h_sub[p * n_sub + q] = h_mo[p * n_basis + q];
    }
  }

  for (int p = 0; p < n_sub; p++) {
    for (int q = 0; q < n_sub; q++) {
      for (int r = 0; r < n_sub; r++) {
        for (int s = 0; s < n_sub; s++) {
          eri_sub[(((size_t)p * n_sub + q) * (size_t)n_sub + r) *
                      (size_t)n_sub +
                  s] = EMO_FULL(p, q, r, s);
        }
      }
    }
  }
#undef EMO_FULL
}

// Builds the active-space Hamiltonian (frozen core folded into an
// effective one-electron operator), diagonalizes it in n_electrons_active
// sector, and constructs the 1-/2-RDM from resulting ground-state CI vector via
// direct fermionic operators. Cross-checked: reconstructing energy as \sum
// h_eff*D1 + 1/2 \sum eri_active*D2 (plus frozen-core/nuclear constant) exactly
// reproduces same ground energy to machine precision
static int casscf_active_space_ci(int n_basis, int n_frozen, int n_active,
                                  int n_electrons_active, const double *h_mo,
                                  const double *eri_mo,
                                  double nuclear_repulsion,
                                  casscf_active_ci_t *out) {
  int n_sub = n_frozen + n_active;

  double *h_sub = malloc((size_t)n_sub * n_sub * sizeof(double));
  double *eri_sub =
      malloc((size_t)n_sub * n_sub * n_sub * n_sub * sizeof(double));
  if (!h_sub || !eri_sub) {
    free(h_sub);
    free(eri_sub);

    return 0;
  }

  casscf_slice_integrals(n_basis, n_sub, h_mo, eri_mo, h_sub, eri_sub);

  cmatrix_t *H = second_quant_build_molecular_hamiltonian_frozen_core(
      n_sub, n_frozen, h_sub, eri_sub, nuclear_repulsion);

  free(h_sub);
  free(eri_sub);
  if (!H) {
    return 0;
  }

  int n_modes = 2 * n_active;
  int full_dim = H->nrows;

  int dim = 0;
  for (int state = 0; state < full_dim; state++) {
    if (__builtin_popcount((unsigned)state) == n_electrons_active) {
      dim++;
    }
  }

  if (dim == 0) {
    cmatrix_free(H);

    return 0;
  }

  int *basis_states = malloc((size_t)dim * sizeof(int));
  int *state_to_local = malloc((size_t)full_dim * sizeof(int));
  if (!basis_states || !state_to_local) {
    free(basis_states);
    free(state_to_local);
    cmatrix_free(H);

    return 0;
  }

  for (int i = 0; i < full_dim; i++) {
    state_to_local[i] = -1;
  }

  int count = 0;
  for (int state = 0; state < full_dim; state++) {
    if (__builtin_popcount((unsigned)state) == n_electrons_active) {
      basis_states[count] = state;
      state_to_local[state] = count;
      count++;
    }
  }

  if (count != dim) {
    free(basis_states);
    free(state_to_local);
    cmatrix_free(H);

    return 0;
  }

  cmatrix_t *H_sector = cmatrix_alloc(dim, dim);
  if (!H_sector) {
    free(basis_states);
    free(state_to_local);
    cmatrix_free(H);

    return 0;
  }

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      CMAT(H_sector, i, j) = CMAT(H, basis_states[i], basis_states[j]);
    }
  }

  cmatrix_free(H);

  eigen_t *eig = cmatrix_eigh_complex(H_sector);
  cmatrix_free(H_sector);
  if (!eig) {
    free(basis_states);
    free(state_to_local);

    return 0;
  }

  double ground_energy = eig->eigenvalues[0];

  double *psi = malloc((size_t)dim * sizeof(double));
  if (!psi) {
    eigen_free(eig);
    free(basis_states);
    free(state_to_local);

    return 0;
  }

  for (int i = 0; i < dim; i++) {
    psi[i] = CMAT(eig->eigenvectors, i, 0).re;
  }

  eigen_free(eig);

  double *D1 = calloc((size_t)n_active * n_active, sizeof(double));
  double *D2 =
      calloc((size_t)n_active * n_active * n_active * n_active, sizeof(double));
  if (!D1 || !D2) {
    free(D1);
    free(D2);
    free(psi);
    free(basis_states);
    free(state_to_local);

    return 0;
  }

  // 1-RDM: D1[P,Q] = <\psi| \sum_\sigma a^+_{P,\sigma} a_{Q,\sigma} |\psi>
  for (int P = 0; P < n_active; P++) {
    for (int Q = 0; Q < n_active; Q++) {
      double v = 0.0;

      for (int sigma = 0; sigma < 2; sigma++) {
        int Ps = 2 * P + sigma, Qs = 2 * Q + sigma;

        for (int j = 0; j < dim; j++) {
          if (fabs(psi[j]) < 1e-14) {
            continue;
          }

          int sign1, sign2;
          int mid =
              casscf_direct_annihilate(basis_states[j], Qs, n_modes, &sign1);
          if (mid < 0) {
            continue;
          }

          int fin = casscf_direct_create(mid, Ps, n_modes, &sign2);
          if (fin < 0) {
            continue;
          }

          int i = state_to_local[fin];
          if (i < 0) {
            continue;
          }

          v += psi[i] * psi[j] * sign1 * sign2;
        }
      }

      D1[P * n_active + Q] = v;
    }
  }

  // 2-RDM: D2[P,Q,R,S] = <\psi| \sum_{\sigma,\tau} a^+_{P,\sigma} a^+_{R,\tau}
  //                            a_{S,\tau} a_{Q,\sigma} |\psi>
  for (int P = 0; P < n_active; P++) {
    for (int Q = 0; Q < n_active; Q++) {
      for (int R = 0; R < n_active; R++) {
        for (int S = 0; S < n_active; S++) {
          double v = 0.0;

          for (int sigma = 0; sigma < 2; sigma++) {
            for (int tau = 0; tau < 2; tau++) {
              int Ps = 2 * P + sigma, Qs = 2 * Q + sigma;
              int Rt = 2 * R + tau, St = 2 * S + tau;

              for (int j = 0; j < dim; j++) {
                if (fabs(psi[j]) < 1e-14) {
                  continue;
                }

                int sign1, sign2, sign3, sign4;
                int m1 = casscf_direct_annihilate(basis_states[j], Qs, n_modes,
                                                  &sign1);
                if (m1 < 0) {
                  continue;
                }

                int m2 = casscf_direct_annihilate(m1, St, n_modes, &sign2);
                if (m2 < 0) {
                  continue;
                }

                int m3 = casscf_direct_create(m2, Rt, n_modes, &sign3);
                if (m3 < 0) {
                  continue;
                }

                int fin = casscf_direct_create(m3, Ps, n_modes, &sign4);
                if (fin < 0) {
                  continue;
                }

                int i = state_to_local[fin];
                if (i < 0) {
                  continue;
                }

                v += psi[i] * psi[j] * sign1 * sign2 * sign3 * sign4;
              }
            }
          }

          D2[((P * n_active + Q) * n_active + R) * n_active + S] = v;
        }
      }
    }
  }

  free(psi);
  free(basis_states);
  free(state_to_local);

  out->energy = ground_energy;
  out->D1 = D1;
  out->D2 = D2;

  return 1;
}

// Extends active-space RDMs to full n_basis-orbital space using MCSCF density
// partition
static void casscf_extend_rdms(int n_basis, int n_frozen, int n_active,
                               const double *D1, const double *D2,
                               double *D_full, double *d_full) {
  for (int i = 0; i < n_basis * n_basis; i++) {
    D_full[i] = 0.0;
  }

  for (int i = 0; i < n_frozen; i++) {
    D_full[i * n_basis + i] = 2.0;
  }

  for (int P = 0; P < n_active; P++) {
    for (int Q = 0; Q < n_active; Q++) {
      D_full[(P + n_frozen) * n_basis + (Q + n_frozen)] = D1[P * n_active + Q];
    }
  }

  size_t n4 = (size_t)n_basis * n_basis * n_basis * n_basis;
  for (size_t i = 0; i < n4; i++) {
    d_full[i] = 0.0;
  }

#define DFULL(p, q, r, s)                                                      \
  d_full[(((size_t)(p) * n_basis + (q)) * (size_t)n_basis + (r)) *             \
             (size_t)n_basis +                                                 \
         (s)]

  for (int i = 0; i < n_frozen; i++) {
    for (int j = 0; j < n_frozen; j++) {
      for (int k = 0; k < n_frozen; k++) {
        for (int l = 0; l < n_frozen; l++) {
          DFULL(i, j, k, l) +=
              4.0 * (i == j) * (k == l) - 2.0 * (i == l) * (j == k);
        }
      }
    }
  }

  for (int i = 0; i < n_frozen; i++) {
    for (int j = 0; j < n_frozen; j++) {
      for (int T = 0; T < n_active; T++) {
        for (int U = 0; U < n_active; U++) {
          int Ta = T + n_frozen, Ua = U + n_frozen;
          double d1tu = D1[T * n_active + U];
          double coulomb = 2.0 * (i == j) * d1tu;
          double exch = -1.0 * (i == j) * d1tu;

          DFULL(i, j, Ta, Ua) += coulomb;
          DFULL(Ta, Ua, i, j) += coulomb;
          DFULL(i, Ta, Ua, j) += exch;
          DFULL(Ta, i, j, Ua) += exch;
        }
      }
    }
  }

  for (int P = 0; P < n_active; P++) {
    for (int Q = 0; Q < n_active; Q++) {
      for (int R = 0; R < n_active; R++) {
        for (int S = 0; S < n_active; S++) {
          DFULL(P + n_frozen, Q + n_frozen, R + n_frozen, S + n_frozen) =
              D2[((P * n_active + Q) * n_active + R) * n_active + S];
        }
      }
    }
  }

#undef DFULL
}

static double casscf_energy_from_rdm(int n_basis, const double *h_mo,
                                     const double *eri_mo,
                                     double nuclear_repulsion,
                                     const double *D_full,
                                     const double *d_full) {
  double E = nuclear_repulsion;

  for (int p = 0; p < n_basis; p++) {
    for (int q = 0; q < n_basis; q++) {
      E += h_mo[p * n_basis + q] * D_full[p * n_basis + q];
    }
  }

  size_t n4 = (size_t)n_basis * n_basis * n_basis * n_basis;
  double two_e = 0.0;
  for (size_t i = 0; i < n4; i++) {
    two_e += eri_mo[i] * d_full[i];
  }

  E += 0.5 * two_e;

  return E;
}

static void casscf_free_active_ci(casscf_active_ci_t *ci) {
  free(ci->D1);
  free(ci->D2);
}

/*
 * C_new = C @ \exp(K),
 * Where
 *  K is built from a flat n_basis x n_basis (row-major) array of real
 * rotation-generator entries.
 */
static cmatrix_t *casscf_rotate(const cmatrix_t *C, const double *kappa_flat,
                                int n) {
  cmatrix_t *K = cmatrix_alloc(n, n);
  if (!K) {
    return NULL;
  }

  for (int i = 0; i < n * n; i++) {
    K->data[i] = c_real(kappa_flat[i]);
  }

  cmatrix_t *U = casscf_matrix_exp(K);
  cmatrix_free(K);
  if (!U) {
    return NULL;
  }

  cmatrix_t *C_new = cmatrix_multiply(C, U);
  cmatrix_free(U);

  return C_new;
}

casscf_result_t *casscf_run(basis_function_t **basis, int n_basis,
                            const molecule_t *mol, const cmatrix_t *C_initial,
                            int n_frozen, int n_active, int n_electrons_active,
                            double conv_tol, int max_iter) {
  if (!basis || n_basis < 1 || !mol || !C_initial || n_frozen < 0 ||
      n_active < 1 || n_frozen + n_active > n_basis || n_electrons_active < 0 ||
      n_electrons_active > 2 * n_active || conv_tol <= 0.0 || max_iter < 1) {
    return NULL;
  }

  cmatrix_t *h_ao = molecular_core_hamiltonian(basis, n_basis, mol);
  double *eri_ao = molecular_eri_tensor(basis, n_basis);
  if (!h_ao || !eri_ao) {
    cmatrix_free(h_ao);
    free(eri_ao);

    return NULL;
  }

  double nuclear_repulsion = molecule_nuclear_repulsion(mol);

  cmatrix_t *C = cmatrix_copy(C_initial);
  size_t n2 = (size_t)n_basis * n_basis;
  size_t n4 = n2 * n2;

  double *h_mo = malloc(n2 * sizeof(double));
  double *eri_mo = malloc(n4 * sizeof(double));
  double *D_full = calloc(n2, sizeof(double));
  double *d_full = calloc(n4, sizeof(double));
  double *grad = malloc(n2 * sizeof(double));
  double *hess_diag = malloc(n2 * sizeof(double));
  double *kappa_flat = malloc(n2 * sizeof(double));

  if (!C || !h_mo || !eri_mo || !D_full || !d_full || !grad || !hess_diag ||
      !kappa_flat) {
    cmatrix_free(h_ao);
    free(eri_ao);
    cmatrix_free(C);
    free(h_mo);
    free(eri_mo);
    free(D_full);
    free(d_full);
    free(grad);
    free(hess_diag);
    free(kappa_flat);

    return NULL;
  }

  const double H_MIN = 0.02;
  const double FD_EPS = 1e-3;

  double E_ci = 0.0;
  double gnorm = 0.0;
  int iter = 0;
  int converged = 0;
  int failed = 0;

  for (iter = 0; iter < max_iter; iter++) {
    molecular_ao_to_mo(h_ao, eri_ao, C, n_basis, h_mo, eri_mo);

    casscf_active_ci_t ci;
    if (!casscf_active_space_ci(n_basis, n_frozen, n_active, n_electrons_active,
                                h_mo, eri_mo, nuclear_repulsion, &ci)) {
      failed = 1;

      break;
    }
    E_ci = ci.energy;
    casscf_extend_rdms(n_basis, n_frozen, n_active, ci.D1, ci.D2, D_full,
                       d_full);
    casscf_free_active_ci(&ci);

    // Orbital gradient + diagonal Hessian via central finite difference,
    // holding D_full/d_full fixed (see casscf.h).
    for (size_t i = 0; i < n2; i++) {
      grad[i] = 0.0;
      hess_diag[i] = 0.0;
    }

    for (int p = 0; p < n_basis; p++) {
      for (int q = p + 1; q < n_basis; q++) {
        for (size_t i = 0; i < n2; i++) {
          kappa_flat[i] = 0.0;
        }
        kappa_flat[p * n_basis + q] = FD_EPS;
        kappa_flat[q * n_basis + p] = -FD_EPS;

        cmatrix_t *Cp = casscf_rotate(C, kappa_flat, n_basis);
        molecular_ao_to_mo(h_ao, eri_ao, Cp, n_basis, h_mo, eri_mo);
        double E_plus = casscf_energy_from_rdm(
            n_basis, h_mo, eri_mo, nuclear_repulsion, D_full, d_full);

        cmatrix_free(Cp);

        kappa_flat[p * n_basis + q] = -FD_EPS;
        kappa_flat[q * n_basis + p] = FD_EPS;

        cmatrix_t *Cm = casscf_rotate(C, kappa_flat, n_basis);
        molecular_ao_to_mo(h_ao, eri_ao, Cm, n_basis, h_mo, eri_mo);
        double E_minus = casscf_energy_from_rdm(
            n_basis, h_mo, eri_mo, nuclear_repulsion, D_full, d_full);

        cmatrix_free(Cm);

        double g = (E_plus - E_minus) / (2.0 * FD_EPS);
        double h = (E_plus - 2.0 * E_ci + E_minus) / (FD_EPS * FD_EPS);

        grad[p * n_basis + q] = g;
        grad[q * n_basis + p] = -g;
        hess_diag[p * n_basis + q] = h;
        hess_diag[q * n_basis + p] = h;
      }
    }

    // restore h_mo/eri_mo to (unrotated) current-orbital values for any code
    // after this point that might need them
    molecular_ao_to_mo(h_ao, eri_ao, C, n_basis, h_mo, eri_mo);

    gnorm = 0.0;
    for (size_t i = 0; i < n2; i++) {
      gnorm += grad[i] * grad[i];
    }

    gnorm = sqrt(gnorm);

    if (gnorm < conv_tol) {
      converged = 1;

      break;
    }

    // Newton-like step, absolute-value-regularized curvature (always valid
    // descent direction)
    for (int p = 0; p < n_basis; p++) {
      for (int q = 0; q < n_basis; q++) {
        double denom = fmax(fabs(hess_diag[p * n_basis + q]), H_MIN);

        kappa_flat[p * n_basis + q] = -grad[p * n_basis + q] / denom;
      }
    }

    // antisymmetrize (already be very close, but keep numerically exact)
    for (int p = 0; p < n_basis; p++) {
      for (int q = p + 1; q < n_basis; q++) {
        double avg =
            0.5 * (kappa_flat[p * n_basis + q] - kappa_flat[q * n_basis + p]);

        kappa_flat[p * n_basis + q] = avg;
        kappa_flat[q * n_basis + p] = -avg;
      }

      kappa_flat[p * n_basis + p] = 0.0;
    }

    double scale = 1.0;
    cmatrix_t *C_new = NULL;
    for (int attempt = 0; attempt < 40; attempt++) {
      double *scaled = malloc(n2 * sizeof(double));
      if (!scaled) {
        break;
      }

      for (size_t i = 0; i < n2; i++) {
        scaled[i] = scale * kappa_flat[i];
      }

      cmatrix_t *C_trial = casscf_rotate(C, scaled, n_basis);
      free(scaled);
      if (!C_trial) {
        break;
      }

      double *h_trial = malloc(n2 * sizeof(double));
      double *eri_trial = malloc(n4 * sizeof(double));
      if (!h_trial || !eri_trial) {
        free(h_trial);
        free(eri_trial);
        cmatrix_free(C_trial);

        break;
      }

      molecular_ao_to_mo(h_ao, eri_ao, C_trial, n_basis, h_trial, eri_trial);

      casscf_active_ci_t ci_trial;
      int ok = casscf_active_space_ci(n_basis, n_frozen, n_active,
                                      n_electrons_active, h_trial, eri_trial,
                                      nuclear_repulsion, &ci_trial);
      free(h_trial);
      free(eri_trial);

      if (ok) {
        double E_trial = ci_trial.energy;

        casscf_free_active_ci(&ci_trial);

        if (E_trial < E_ci + 1e-10) {
          C_new = C_trial;

          break;
        }
      }

      cmatrix_free(C_trial);
      scale *= 0.5;
    }

    if (!C_new) {
      failed = 1;

      break;
    }

    cmatrix_free(C);
    C = C_new;
  }

  cmatrix_free(h_ao);
  free(eri_ao);
  free(h_mo);
  free(eri_mo);
  free(D_full);
  free(d_full);
  free(grad);
  free(hess_diag);
  free(kappa_flat);

  if (failed) {
    cmatrix_free(C);

    return NULL;
  }

  casscf_result_t *result = malloc(sizeof(casscf_result_t));
  if (!result) {
    cmatrix_free(C);

    return NULL;
  }

  result->total_energy = E_ci;
  result->n_spatial = n_basis;
  result->n_frozen = n_frozen;
  result->n_active = n_active;
  result->n_electrons_active = n_electrons_active;
  result->C = C;
  result->iterations = iter;
  result->converged = converged;
  result->grad_norm = gnorm;

  return result;
}

void casscf_result_free(casscf_result_t *res) {
  if (!res) {
    return;
  }

  cmatrix_free(res->C);
  free(res);
}
