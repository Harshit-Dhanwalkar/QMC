#include "molecular_integrals.h"
#include "../core/complex.h"
#include "../core/matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------
 * Boys function
 * ------------------------------------------------------------------- */

static double boys_series(int n, double x, int nterms) {
  /* F_n(x) = \exp(-x) * \sum_{k = 0}^\inf [(2n-1)!! / (2n+2k+1)!!] (2x)^k
   * NOTE: All terms in the sum are positive, factoring out \exp(-x) makes this stable
   * for large x too. */
  double s = 1.0 / (2 * n + 1);
  double term = s;

  for (int k = 1; k < nterms; k++) {
    term *= (2.0 * x) / (2 * n + 2 * k + 1);
    s += term;

    if (fabs(term) < 1e-17 * fabs(s)) {
      break;
    }
  }

  return exp(-x) * s;
}

void boys_function_array(int nmax, double x, double *F) {
  if (nmax < 0) {
    return;
  }

  F[nmax] = boys_series(nmax, x, 300);
  double ex = exp(-x);

  for (int n = nmax; n > 0; n--) {
    F[n - 1] = (2.0 * x * F[n] + ex) / (2 * n - 1);
  }
}

double boys_function(int n, double x) {
  double *F = malloc((size_t)(n + 1) * sizeof(double));
  if (!F) {
    return 0.0;
  }

  boys_function_array(n, x, F);
  double result = F[n];

  free(F);

  return result;
}

/* ---------------------------------------------------------------------
 * McMurchie-Davidson 1D Hermite expansion coefficients E^{ij}_t
 * ------------------------------------------------------------------- */

/* NOTE: Recursive, matching the validated Python 1:1. i,j are bounded by basis
 * angular momentum (small: s/p/d/f in practice), so unmemoized recursion is not
 * a performance concern here.
 */
static double md_E(int i, int j, int t, double Qx, double a, double b) {
  double p = a + b;
  double q = a * b / p;

  if (t < 0 || t > i + j) {
    return 0.0;
  }

  if (i == 0 && j == 0 && t == 0) {
    return exp(-q * Qx * Qx);
  }

  if (j == 0) {
    double term = 0.0;
    if (i - 1 >= 0) {
      term = (1.0 / (2.0 * p)) * md_E(i - 1, j, t - 1, Qx, a, b) -
             (q * Qx / a) * md_E(i - 1, j, t, Qx, a, b) +
             (t + 1) * md_E(i - 1, j, t + 1, Qx, a, b);
    }

    return term;
  } else {
    double term = (1.0 / (2.0 * p)) * md_E(i, j - 1, t - 1, Qx, a, b) +
                  (q * Qx / b) * md_E(i, j - 1, t, Qx, a, b) +
                  (t + 1) * md_E(i, j - 1, t + 1, Qx, a, b);

    return term;
  }
}

static double overlap_1d(int i, int j, double Qx, double a, double b) {
  double p = a + b;

  return md_E(i, j, 0, Qx, a, b) * sqrt(M_PI / p);
}

static double kinetic_1d(int i, int j, double Qx, double a, double b) {
  double term1 = b * (2 * j + 1) * overlap_1d(i, j, Qx, a, b);
  double term2 = -2.0 * b * b * overlap_1d(i, j + 2, Qx, a, b);
  double term3 =
      (j >= 2) ? -0.5 * j * (j - 1) * overlap_1d(i, j - 2, Qx, a, b) : 0.0;

  return term1 + term2 + term3;
}

/* NOTE: Hermite Coulomb integral R^n_{tuv}(p, PC), recursive per Helgaker et
 * al. Fvals must already hold Boys ladder F[0..Lmax] for this center pair at
 * argument p*|PC|^2 (computed once per primitive pair, not once per (t,u,v)).
 */
static double md_R(int t, int u, int v, int n, double p, double PCx, double PCy,
                   double PCz, double *Fvals) {
  if (t == 0 && u == 0 && v == 0) {
    return pow(-2.0 * p, n) * Fvals[n];
  }

  if (t > 0) {
    double val = 0.0;

    if (t > 1) {
      val += (t - 1) * md_R(t - 2, u, v, n + 1, p, PCx, PCy, PCz, Fvals);
    }
    val += PCx * md_R(t - 1, u, v, n + 1, p, PCx, PCy, PCz, Fvals);

    return val;
  }

  if (u > 0) {
    double val = 0.0;

    if (u > 1) {
      val += (u - 1) * md_R(t, u - 2, v, n + 1, p, PCx, PCy, PCz, Fvals);
    }
    val += PCy * md_R(t, u - 1, v, n + 1, p, PCx, PCy, PCz, Fvals);

    return val;
  }

  {
    double val = 0.0;

    if (v > 1) {
      val += (v - 1) * md_R(t, u, v - 2, n + 1, p, PCx, PCy, PCz, Fvals);
    }
    val += PCz * md_R(t, u, v - 1, n + 1, p, PCx, PCy, PCz, Fvals);

    return val;
  }
}

/* ---------------------------------------------------------------------
 * Primitive-pair integrals (single Gaussian x single Gaussian)
 * ------------------------------------------------------------------- */

static double prim_overlap(const gto_primitive_t *A, const double Ac[3],
                           const gto_primitive_t *B, const double Bc[3]) {
  double a = A->alpha, b = B->alpha;
  double Sx = overlap_1d(A->l, B->l, Ac[0] - Bc[0], a, b);
  double Sy = overlap_1d(A->m, B->m, Ac[1] - Bc[1], a, b);
  double Sz = overlap_1d(A->n, B->n, Ac[2] - Bc[2], a, b);

  return Sx * Sy * Sz;
}

static double prim_kinetic(const gto_primitive_t *A, const double Ac[3],
                           const gto_primitive_t *B, const double Bc[3]) {
  double a = A->alpha, b = B->alpha;
  double Qx = Ac[0] - Bc[0], Qy = Ac[1] - Bc[1], Qz = Ac[2] - Bc[2];
  double Sx = overlap_1d(A->l, B->l, Qx, a, b);
  double Sy = overlap_1d(A->m, B->m, Qy, a, b);
  double Sz = overlap_1d(A->n, B->n, Qz, a, b);
  double Tx = kinetic_1d(A->l, B->l, Qx, a, b);
  double Ty = kinetic_1d(A->m, B->m, Qy, a, b);
  double Tz = kinetic_1d(A->n, B->n, Qz, a, b);

  return Tx * Sy * Sz + Sx * Ty * Sz + Sx * Sy * Tz;
}

static double prim_nuclear(const gto_primitive_t *A, const double Ac[3],
                           const gto_primitive_t *B, const double Bc[3],
                           const double C[3]) {
  double a = A->alpha, b = B->alpha, p = a + b;
  double P[3] = {(a * Ac[0] + b * Bc[0]) / p, (a * Ac[1] + b * Bc[1]) / p,
                 (a * Ac[2] + b * Bc[2]) / p};
  double PC[3] = {P[0] - C[0], P[1] - C[1], P[2] - C[2]};
  double RPC2 = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

  int Lmax = A->l + B->l + A->m + B->m + A->n + B->n;
  double *Fvals = malloc((size_t)(Lmax + 1) * sizeof(double));
  if (!Fvals) {
    return 0.0;
  }
  boys_function_array(Lmax, p * RPC2, Fvals);

  double total = 0.0;
  for (int t = 0; t <= A->l + B->l; t++) {
    double Ex = md_E(A->l, B->l, t, Ac[0] - Bc[0], a, b);
    if (Ex == 0.0) {
      continue;
    }

    for (int u = 0; u <= A->m + B->m; u++) {
      double Ey = md_E(A->m, B->m, u, Ac[1] - Bc[1], a, b);
      if (Ey == 0.0) {
        continue;
      }

      for (int v = 0; v <= A->n + B->n; v++) {
        double Ez = md_E(A->n, B->n, v, Ac[2] - Bc[2], a, b);
        if (Ez == 0.0) {
          continue;
        }

        total += Ex * Ey * Ez * md_R(t, u, v, 0, p, PC[0], PC[1], PC[2], Fvals);
      }
    }
  }

  free(Fvals);

  return total * 2.0 * M_PI / p;
}

static double prim_eri(const gto_primitive_t *A, const double Ac[3],
                       const gto_primitive_t *B, const double Bc[3],
                       const gto_primitive_t *Cp, const double Cc[3],
                       const gto_primitive_t *D, const double Dc[3]) {
  double a = A->alpha, b = B->alpha, c = Cp->alpha, d = D->alpha;
  double p = a + b, q = c + d;
  double P[3] = {(a * Ac[0] + b * Bc[0]) / p, (a * Ac[1] + b * Bc[1]) / p,
                 (a * Ac[2] + b * Bc[2]) / p};
  double Q[3] = {(c * Cc[0] + d * Dc[0]) / q, (c * Cc[1] + d * Dc[1]) / q,
                 (c * Cc[2] + d * Dc[2]) / q};
  double alpha = p * q / (p + q);
  double PQ[3] = {P[0] - Q[0], P[1] - Q[1], P[2] - Q[2]};
  double RPQ2 = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

  int Lmax = (A->l + B->l + A->m + B->m + A->n + B->n) +
             (Cp->l + D->l + Cp->m + D->m + Cp->n + D->n);
  double *Fvals = malloc((size_t)(Lmax + 1) * sizeof(double));
  if (!Fvals) {
    return 0.0;
  }
  boys_function_array(Lmax, alpha * RPQ2, Fvals);

  double total = 0.0;
  for (int t1 = 0; t1 <= A->l + B->l; t1++) {
    double Ex1 = md_E(A->l, B->l, t1, Ac[0] - Bc[0], a, b);

    if (Ex1 == 0.0) {
      continue;
    }

    for (int u1 = 0; u1 <= A->m + B->m; u1++) {
      double Ey1 = md_E(A->m, B->m, u1, Ac[1] - Bc[1], a, b);
      if (Ey1 == 0.0) {
        continue;
      }

      for (int v1 = 0; v1 <= A->n + B->n; v1++) {
        double Ez1 = md_E(A->n, B->n, v1, Ac[2] - Bc[2], a, b);
        if (Ez1 == 0.0) {
          continue;
        }

        double E1 = Ex1 * Ey1 * Ez1;

        for (int t2 = 0; t2 <= Cp->l + D->l; t2++) {
          double Ex2 = md_E(Cp->l, D->l, t2, Cc[0] - Dc[0], c, d);
          if (Ex2 == 0.0) {
            continue;
          }

          for (int u2 = 0; u2 <= Cp->m + D->m; u2++) {
            double Ey2 = md_E(Cp->m, D->m, u2, Cc[1] - Dc[1], c, d);
            if (Ey2 == 0.0) {
              continue;
            }

            for (int v2 = 0; v2 <= Cp->n + D->n; v2++) {
              double Ez2 = md_E(Cp->n, D->n, v2, Cc[2] - Dc[2], c, d);
              if (Ez2 == 0.0) {
                continue;
              }

              double E2 = Ex2 * Ey2 * Ez2;

              double sign = ((t2 + u2 + v2) % 2 == 0) ? 1.0 : -1.0;
              total += E1 * E2 * sign *
                       md_R(t1 + t2, u1 + u2, v1 + v2, 0, alpha, PQ[0], PQ[1],
                            PQ[2], Fvals);
            }
          }
        }
      }
    }
  }

  free(Fvals);

  return total * 2.0 * pow(M_PI, 2.5) / (p * q * sqrt(p + q));
}

/* ---------------------------------------------------------------------
 * Contracted (basis_function_t) integrals: sum over primitive pairs
 * ------------------------------------------------------------------- */

double gto_overlap(const basis_function_t *A, const basis_function_t *B) {
  double s = 0.0;
  for (int i = 0; i < A->n_primitives; i++) {
    gto_primitive_t pa = {A->l, A->m, A->n, A->exponents[i]};

    for (int j = 0; j < B->n_primitives; j++) {
      gto_primitive_t pb = {B->l, B->m, B->n, B->exponents[j]};
      s += A->coefficients[i] * B->coefficients[j] *
           prim_overlap(&pa, A->center, &pb, B->center);
    }
  }

  return s;
}

double gto_kinetic(const basis_function_t *A, const basis_function_t *B) {
  double s = 0.0;
  for (int i = 0; i < A->n_primitives; i++) {
    gto_primitive_t pa = {A->l, A->m, A->n, A->exponents[i]};

    for (int j = 0; j < B->n_primitives; j++) {
      gto_primitive_t pb = {B->l, B->m, B->n, B->exponents[j]};
      s += A->coefficients[i] * B->coefficients[j] *
           prim_kinetic(&pa, A->center, &pb, B->center);
    }
  }

  return s;
}

double gto_nuclear_attraction(const basis_function_t *A,
                              const basis_function_t *B,
                              const double center[3]) {
  double s = 0.0;
  for (int i = 0; i < A->n_primitives; i++) {
    gto_primitive_t pa = {A->l, A->m, A->n, A->exponents[i]};

    for (int j = 0; j < B->n_primitives; j++) {
      gto_primitive_t pb = {B->l, B->m, B->n, B->exponents[j]};
      s += A->coefficients[i] * B->coefficients[j] *
           prim_nuclear(&pa, A->center, &pb, B->center, center);
    }
  }

  return s;
}

double gto_eri(const basis_function_t *A, const basis_function_t *B,
               const basis_function_t *C, const basis_function_t *D) {
  double s = 0.0;
  for (int i = 0; i < A->n_primitives; i++) {
    gto_primitive_t pa = {A->l, A->m, A->n, A->exponents[i]};
    for (int j = 0; j < B->n_primitives; j++) {
      gto_primitive_t pb = {B->l, B->m, B->n, B->exponents[j]};
      for (int k = 0; k < C->n_primitives; k++) {
        gto_primitive_t pc = {C->l, C->m, C->n, C->exponents[k]};
        for (int m = 0; m < D->n_primitives; m++) {
          gto_primitive_t pd = {D->l, D->m, D->n, D->exponents[m]};
          s += A->coefficients[i] * B->coefficients[j] * C->coefficients[k] *
               D->coefficients[m] *
               prim_eri(&pa, A->center, &pb, B->center, &pc, C->center, &pd,
                        D->center);
        }
      }
    }
  }

  return s;
}

/* ---------------------------------------------------------------------
 * basis_function_t / molecule_t construction and normalization
 * ------------------------------------------------------------------- */

basis_function_t *basis_function_alloc(int l, int m, int n,
                                       const double center[3], int n_primitives,
                                       const double *exponents,
                                       const double *raw_coefficients) {
  if (n_primitives <= 0) {
    return NULL;
  }

  basis_function_t *bf = malloc(sizeof(basis_function_t));
  if (!bf) {
    return NULL;
  }

  bf->l = l;
  bf->m = m;
  bf->n = n;
  bf->center[0] = center[0];
  bf->center[1] = center[1];
  bf->center[2] = center[2];
  bf->n_primitives = n_primitives;
  bf->exponents = malloc((size_t)n_primitives * sizeof(double));
  bf->coefficients = malloc((size_t)n_primitives * sizeof(double));
  if (!bf->exponents || !bf->coefficients) {
    free(bf->exponents);
    free(bf->coefficients);
    free(bf);

    return NULL;
  }

  memcpy(bf->exponents, exponents, (size_t)n_primitives * sizeof(double));
  memcpy(bf->coefficients, raw_coefficients,
         (size_t)n_primitives * sizeof(double));

  return bf;
}

void basis_function_free(basis_function_t *bf) {
  if (!bf) {
    return;
  }

  free(bf->exponents);
  free(bf->coefficients);
  free(bf);
}

void molint_normalize_contraction(basis_function_t *bf) {
  /* Primitive Cartesian-Gaussian normalization constant for :
   *    x^l y^m z^n \exp(-\alpha * r^2), general (l,m,n)
   *  NOTE: (Helgaker et al. eq. 6.24-ish; for s functions this reduces to :
   *    (2 * a / \pi)^{3 / 4})
   */
  int L = bf->l + bf->m + bf->n;
  for (int i = 0; i < bf->n_primitives; i++) {
    double a = bf->exponents[i];
    double num = pow(2.0 * a / M_PI, 1.5) * pow(4.0 * a, L);
    double denom = 1.0;

    // (2l-1)!!(2m-1)!!(2n-1)!!
    for (int k = bf->l; k > 1; k -= 2) {
      denom *= k;
    }

    for (int k = bf->m; k > 1; k -= 2) {
      denom *= k;
    }

    for (int k = bf->n; k > 1; k -= 2) {
      denom *= k;
    }

    bf->coefficients[i] *= sqrt(num / denom);
  }

  // Rescale whole contraction so <bf|bf> = 1, using gto_overlap
  double self_overlap = gto_overlap(bf, bf);
  if (self_overlap > 0.0) {
    double norm = 1.0 / sqrt(self_overlap);

    for (int i = 0; i < bf->n_primitives; i++) {
      bf->coefficients[i] *= norm;
    }
  }
}

molecule_t *molecule_alloc(int n_atoms, const double *charge,
                           const double center[][3]) {
  if (n_atoms <= 0) {
    return NULL;
  }

  molecule_t *mol = malloc(sizeof(molecule_t));
  if (!mol) {
    return NULL;
  }

  mol->n_atoms = n_atoms;
  mol->charge = malloc((size_t)n_atoms * sizeof(double));
  mol->center = malloc((size_t)n_atoms * sizeof(*mol->center));
  if (!mol->charge || !mol->center) {
    free(mol->charge);
    free(mol->center);
    free(mol);

    return NULL;
  }

  memcpy(mol->charge, charge, (size_t)n_atoms * sizeof(double));
  memcpy(mol->center, center, (size_t)n_atoms * sizeof(*mol->center));

  return mol;
}

void molecule_free(molecule_t *mol) {
  if (!mol) {
    return;
  }

  free(mol->charge);
  free(mol->center);
  free(mol);
}

double molecule_nuclear_repulsion(const molecule_t *mol) {
  double e = 0.0;
  for (int i = 0; i < mol->n_atoms; i++) {
    for (int j = i + 1; j < mol->n_atoms; j++) {
      double dx = mol->center[i][0] - mol->center[j][0];
      double dy = mol->center[i][1] - mol->center[j][1];
      double dz = mol->center[i][2] - mol->center[j][2];
      double r = sqrt(dx * dx + dy * dy + dz * dz);

      if (r > 0.0) {
        e += mol->charge[i] * mol->charge[j] / r;
      }
    }
  }

  return e;
}

/* ---------------------------------------------------------------------
 * Whole-basis builders
 * ------------------------------------------------------------------- */

cmatrix_t *molecular_overlap_matrix(basis_function_t **basis, int n_basis) {
  if (n_basis <= 0) {
    return NULL;
  }

  cmatrix_t *S = cmatrix_alloc(n_basis, n_basis);
  if (!S) {
    return NULL;
  }

  for (int i = 0; i < n_basis; i++) {
    for (int j = 0; j < n_basis; j++) {
      double v = gto_overlap(basis[i], basis[j]);

      CMAT(S, i, j) = (complex_t){v, 0.0};
    }
  }

  return S;
}

cmatrix_t *molecular_kinetic_matrix(basis_function_t **basis, int n_basis) {
  if (n_basis <= 0) {
    return NULL;
  }

  cmatrix_t *T = cmatrix_alloc(n_basis, n_basis);
  if (!T) {
    return NULL;
  }

  for (int i = 0; i < n_basis; i++) {
    for (int j = 0; j < n_basis; j++) {
      double v = gto_kinetic(basis[i], basis[j]);

      CMAT(T, i, j) = (complex_t){v, 0.0};
    }
  }

  return T;
}

cmatrix_t *molecular_core_hamiltonian(basis_function_t **basis, int n_basis,
                                      const molecule_t *mol) {
  if (n_basis <= 0) {
    return NULL;
  }

  cmatrix_t *H = cmatrix_alloc(n_basis, n_basis);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < n_basis; i++) {
    for (int j = 0; j < n_basis; j++) {
      double v = gto_kinetic(basis[i], basis[j]);

      for (int A = 0; A < mol->n_atoms; A++) {
        v -= mol->charge[A] *
             gto_nuclear_attraction(basis[i], basis[j], mol->center[A]);
      }

      CMAT(H, i, j) = (complex_t){v, 0.0};
    }
  }

  return H;
}

double *molecular_eri_tensor(basis_function_t **basis, int n_basis) {
  if (n_basis <= 0) {
    return NULL;
  }

  size_t n4 = (size_t)n_basis * n_basis * n_basis * n_basis;
  double *eri = malloc(n4 * sizeof(double));
  if (!eri) {
    return NULL;
  }

  for (size_t k = 0; k < n4; k++) {
    eri[k] = NAN; // not yet computed
  }

  for (int i = 0; i < n_basis; i++) {
    for (int j = 0; j <= i; j++) {
      for (int k = 0; k < n_basis; k++) {
        for (int l = 0; l <= k; l++) {
          /* (ij|kl), computing each symmetry-distinct integral once via
           * canonical (i>=j, k>=l, (ij)>=(kl) as pair-index) ordering, then
           * copying to all 8 equivalent index permutations. */
          size_t ij = (size_t)i * (i + 1) / 2 + j;
          size_t kl = (size_t)k * (k + 1) / 2 + l;
          if (kl > ij) {
            continue;
          }

          double val = gto_eri(basis[i], basis[j], basis[k], basis[l]);
          int idx_i[8] = {i, j, i, j, k, l, k, l};
          int idx_j[8] = {j, i, j, i, l, k, l, k};
          int idx_k[8] = {k, k, l, l, i, i, j, j};
          int idx_l[8] = {l, l, k, k, j, j, i, i};
          for (int s = 0; s < 8; s++) {
            MOLINT_ERI(eri, n_basis, idx_i[s], idx_j[s], idx_k[s], idx_l[s]) =
                val;
          }
        }
      }
    }
  }

  return eri;
}

/* ---------------------------------------------------------------------
 * STO-3G hydrogen basis builder
 * ------------------------------------------------------------------- */

basis_function_t *molint_basis_sto3g_h(const double center[3]) {
  /* NOTE: Published STO-3G H 1s parameters (Szabo & Ostlund Table 3.7 /
   * EMSL/Gaussian94 STO-3G basis set for hydrogen). */
  static const double alphas[3] = {3.42525091, 0.62391373, 0.16885540};
  static const double raw_coeffs[3] = {0.15432897, 0.53532814, 0.44463454};

  basis_function_t *bf =
      basis_function_alloc(0, 0, 0, center, 3, alphas, raw_coeffs);
  if (!bf) {
    return NULL;
  }

  molint_normalize_contraction(bf);

  return bf;
}
