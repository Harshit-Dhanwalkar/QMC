/*
 * UCC ansatz implementation (dense small-system reference)
 *
 * Generator construction (ucc_build_generator):
 *   For each single (i, a), term = a_a^+ a_i and term_dag = a_i^+ a_a are
 *   built from JW matrices (jw_creation_operator / jw_annihilation_operator)
 *   Contribution to kappa is \theta_s * (term - term_dag), which is
 *   anti-Hermitian by construction. For each double (i, j, a, b), term
 *   is a_a^+ a_b^+ a_j a_i and its adjoint is a_i^+ a_j^+ a_b a_a; again
 *   \theta_d * (term - term_dag) is anti-Hermitian. \kappa is a single
 *   accumulated dense matrix; intermediate term/term_dag/tensor-product
 *   matrices are freed as soon as they've been folded in
 *   Index validation happens before any operator is built, so an invalid (i, a)
 *   or (i, j, a, b) always returns NULL with no partial state
 *
 * Exponential (ucc_prepare_state):
 *   H = i * kappa is Hermitian when \kappa is anti-Hermitian, so we diagonalize
 *   H instead of kappa - this lets us reuse cmatrix_eigh_complex directly
 *   (its input contract is Hermitian) and gives real eigenvalues \lambda_k
 *   Then \exp(\kappa) = \exp(-i*H) = V * diag(\exp(-i * \lambda_k)) * V^\dagger
 *   We apply this to |reference> as:
 *       coeffs = V^\dagger @ reference      (project)
 *       coeffs[k] *= \exp(-i * \lambda_k)   (phase in eigenbasis)
 *       \psi = V @ coeffs                   (back-transform)
 *   which is O(dim^2) once eigenvectors are known, instead of O(dim^3) to form
 *   \exp(\kappa) explicitly. Result is normalized with cvector_normalize
 *   afterward: \exp(\kappa) is exactly unitary only in exact arithmetic;
 *   rounding in eigendecomposition can leave norm at 1 \pm 1e-15 and
 *   re-normalizing keeps downstream overlap computations consistent with each
 *   other
 *
 * CCSD amplitude conversion (ucc_excitations_from_ccsd_amplitudes):
 *   t1[i, a] and t2[i, j, a, b] are already "converged amplitude" for single (i
 *   -> a) and double (i,j -> a,b) excitations. Generator uses these directly as
 *   \theta_s / \theta_d, which is standard "UCC from CC amplitudes"
 *   initialization; it is not variational optimum of UCC, just a good starting
 *   point. Occupied and virtual indices are taken from amp->occ[]/amp->virt[]
 *   (they are not necessarily contiguous 0..no-1 / no..nso-1) so returned
 *   indices are always valid spin-orbital labels even for restricted
 *   calculations where CCSD is run in a symmetry-adapted (frozen-core) space
 *
 *   Doubles enumeration uses i < j (occupied) and a < b (virtual) to avoid
 *   double-counting antisymmetric pair; single t2 entry t2[i, j, a, b] carries
 *   full amplitude, so (j, i, b, a) variant would just be its negation and must
 *   not be added separately
 *
 *   All four output pointers are set on success. If n_doubles == 0 caller gets
 *   NULL doubles/\theta_d with *n_doubles_out = 0 - this is a legitimate case
 *   (e.g. one virtual orbital) and is handled uniformly by ucc_build_generator,
 *   which accepts NULL list when count == 0
 */
#include "ucc.h"
#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../core/vector.h"
#include "ccsd.h"
#include "second_quant.h"
#include <math.h>
#include <stdlib.h>

#define IDX2(arr, nso, p, q) (arr)[(size_t)(p) * (nso) + (q)]
#define IDX4(arr, nso, p, q, r, s)                                             \
  (arr)[(((size_t)(p) * (nso) + (q)) * (nso) + (r)) * (nso) + (s)]

cmatrix_t *ucc_build_generator(int n_modes, const ucc_single_t *singles,
                               const double *theta_s, int n_singles,
                               const ucc_double_t *doubles,
                               const double *theta_d, int n_doubles) {
  if (n_modes < 1 || (n_singles > 0 && (!singles || !theta_s)) ||
      (n_doubles > 0 && (!doubles || !theta_d))) {
    return NULL;
  }

  int dim = 1 << n_modes;
  cmatrix_t *kappa = cmatrix_alloc(dim, dim);
  if (!kappa) {
    return NULL;
  }

  for (int i = 0; i < dim * dim; i++) {
    kappa->data[i] = c_zero();
  }

  for (int k = 0; k < n_singles; k++) {
    int i = singles[k].i, a = singles[k].a;
    if (i < 0 || i >= n_modes || a < 0 || a >= n_modes) {
      cmatrix_free(kappa);

      return NULL;
    }

    cmatrix_t *ai_dag = jw_creation_operator(a, n_modes);
    cmatrix_t *ii_ann = jw_annihilation_operator(i, n_modes);
    cmatrix_t *term = cmatrix_multiply(ai_dag, ii_ann); // a_a^+ a_i

    cmatrix_free(ai_dag);
    cmatrix_free(ii_ann);

    cmatrix_t *term_dag = cmatrix_adjoint(term); // a_i^+ a_a

    for (int e = 0; e < dim * dim; e++) {
      complex_t contrib =
          c_scale(c_sub(term->data[e], term_dag->data[e]), theta_s[k]);
      kappa->data[e] = c_add(kappa->data[e], contrib);
    }

    cmatrix_free(term);
    cmatrix_free(term_dag);
  }

  for (int k = 0; k < n_doubles; k++) {
    int i = doubles[k].i, j = doubles[k].j, a = doubles[k].a, b = doubles[k].b;
    if (i < 0 || i >= n_modes || j < 0 || j >= n_modes || a < 0 ||
        a >= n_modes || b < 0 || b >= n_modes) {
      cmatrix_free(kappa);

      return NULL;
    }

    cmatrix_t *ad = jw_creation_operator(a, n_modes);
    cmatrix_t *bd = jw_creation_operator(b, n_modes);
    cmatrix_t *aj = jw_annihilation_operator(j, n_modes);
    cmatrix_t *ai = jw_annihilation_operator(i, n_modes);

    cmatrix_t *ad_bd = cmatrix_multiply(ad, bd);
    cmatrix_t *aj_ai = cmatrix_multiply(aj, ai);
    cmatrix_t *term = cmatrix_multiply(ad_bd, aj_ai); // a_a^+ a_b^+ a_j a_i

    cmatrix_free(ad);
    cmatrix_free(bd);
    cmatrix_free(aj);
    cmatrix_free(ai);
    cmatrix_free(ad_bd);
    cmatrix_free(aj_ai);

    cmatrix_t *term_dag = cmatrix_adjoint(term); // a_i^+ a_j^+ a_b a_a

    for (int e = 0; e < dim * dim; e++) {
      complex_t contrib =
          c_scale(c_sub(term->data[e], term_dag->data[e]), theta_d[k]);

      kappa->data[e] = c_add(kappa->data[e], contrib);
    }

    cmatrix_free(term);
    cmatrix_free(term_dag);
  }

  return kappa;
}

cvector_t *ucc_prepare_state(const cmatrix_t *generator,
                             const cvector_t *reference) {
  if (!generator || !reference || generator->nrows != generator->ncols ||
      generator->nrows != reference->n) {
    return NULL;
  }

  int dim = generator->nrows;

  // H = i * generator, Hermitian since generator is anti-Hermitian
  cmatrix_t *H = cmatrix_alloc(dim, dim);
  if (!H) {
    return NULL;
  }

  for (int e = 0; e < dim * dim; e++) {
    H->data[e] = c_mul(c_new(0.0, 1.0), generator->data[e]);
  }

  eigen_t *eig = cmatrix_eigh_complex(H);
  cmatrix_free(H);
  if (!eig) {
    return NULL;
  }

  // coeffs = V^\dagger @ reference
  cvector_t *coeffs = cvector_alloc(dim);
  if (!coeffs) {
    eigen_free(eig);

    return NULL;
  }

  for (int k = 0; k < dim; k++) {
    complex_t s = c_zero();

    for (int p = 0; p < dim; p++) {
      s = c_add(
          s, c_mul(c_conj(CMAT(eig->eigenvectors, p, k)), reference->data[p]));
    }

    // \exp(-i * \lambda_k) applied in eigenbasis
    complex_t phase =
        c_new(cos(-eig->eigenvalues[k]), sin(-eig->eigenvalues[k]));

    coeffs->data[k] = c_mul(s, phase);
  }

  // psi = V @ coeffs
  cvector_t *psi = cvector_alloc(dim);
  if (!psi) {
    cvector_free(coeffs);
    eigen_free(eig);

    return NULL;
  }

  for (int p = 0; p < dim; p++) {
    complex_t s = c_zero();

    for (int k = 0; k < dim; k++) {
      s = c_add(s, c_mul(CMAT(eig->eigenvectors, p, k), coeffs->data[k]));
    }

    psi->data[p] = s;
  }

  cvector_free(coeffs);
  eigen_free(eig);
  cvector_normalize(psi);

  return psi;
}

int ucc_excitations_from_ccsd_amplitudes(
    const ccsd_amplitudes_t *amp, ucc_single_t **singles_out,
    double **theta_s_out, int *n_singles_out, ucc_double_t **doubles_out,
    double **theta_d_out, int *n_doubles_out) {
  if (!amp || !singles_out || !theta_s_out || !n_singles_out || !doubles_out ||
      !theta_d_out || !n_doubles_out) {
    return 0;
  }

  int nso = amp->nso;
  int no = amp->nocc, nv = amp->nvirt;
  const int *occ = amp->occ, *virt = amp->virt;

  int n_singles = no * nv;
  ucc_single_t *singles = malloc((size_t)n_singles * sizeof(ucc_single_t));
  double *theta_s = malloc((size_t)n_singles * sizeof(double));
  if (!singles || !theta_s) {
    free(singles);
    free(theta_s);

    return 0;
  }

  int k = 0;
  for (int ii = 0; ii < no; ii++) {
    for (int ai = 0; ai < nv; ai++) {
      singles[k].i = occ[ii];
      singles[k].a = virt[ai];
      theta_s[k] = IDX2(amp->t1, nso, occ[ii], virt[ai]);
      k++;
    }
  }

  int n_pairs_occ = no * (no - 1) / 2;
  int n_pairs_virt = nv * (nv - 1) / 2;
  int n_doubles = n_pairs_occ * n_pairs_virt;

  ucc_double_t *doubles = NULL;
  double *theta_d = NULL;
  if (n_doubles > 0) {
    doubles = malloc((size_t)n_doubles * sizeof(ucc_double_t));
    theta_d = malloc((size_t)n_doubles * sizeof(double));

    if (!doubles || !theta_d) {
      free(singles);
      free(theta_s);
      free(doubles);
      free(theta_d);

      return 0;
    }
  }

  int kd = 0;
  for (int ii = 0; ii < no; ii++) {
    for (int ji = ii + 1; ji < no; ji++) {
      for (int ai = 0; ai < nv; ai++) {
        for (int bi = ai + 1; bi < nv; bi++) {
          int i = occ[ii];
          int j = occ[ji];
          int a = virt[ai];
          int b = virt[bi];

          doubles[kd].i = i;
          doubles[kd].j = j;
          doubles[kd].a = a;
          doubles[kd].b = b;

          theta_d[kd] = IDX4(amp->t2, nso, i, j, a, b);
          kd++;
        }
      }
    }
  }

  *singles_out = singles;
  *theta_s_out = theta_s;
  *n_singles_out = n_singles;
  *doubles_out = doubles;
  *theta_d_out = theta_d;
  *n_doubles_out = n_doubles;

  return 1;
}
