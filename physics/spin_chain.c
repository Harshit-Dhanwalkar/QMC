#include "spin_chain.h"
#include "../core/complex.h"
#include "../core/sparse.h"
#include "../core/vector.h"
#include "/home/harshitpd/Documents/GITHUB/QMC/core/matrix.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* bit-level helpers
 *
 * NOTE: Bit j = 0 -> spin up (S^z_j = +1/2); bit j = 1 -> spin down (S^z_j =
 * -1/2). translate() implements one step of the cyclic shift T: site j -> site
 * j+1 (mod N), i.e. it moves every bit one position to the left, wrapping
 * the bit that falls off the top back into bit 0.
 */

static inline uint64_t translate(uint64_t s, int N, uint64_t mask) {
  uint64_t top = (s >> (N - 1)) & 1ULL;

  return ((s << 1) | top) & mask;
}

static inline double sz_bit(uint64_t s, int j) {
  return ((s >> j) & 1ULL) ? -0.5 : 0.5;
}

/* Find the representative (numerically smallest bitstring in the
 * translation orbit) of s, and l such that translate^l(rep) == s.
 *
 * This scans the whole N-element orbit, so sector construction and
 * Hamiltonian assembly
 * WARN: sector dimension itself is bottleneck, not this scan.
 */
static uint64_t find_representative(uint64_t s, int N, uint64_t mask,
                                    int *l_out) {
  uint64_t best = s;
  int lbest = 0;
  uint64_t t = s;

  for (int l = 1; l < N; l++) {
    t = translate(t, N, mask);
    if (t < best) {
      best = t;
      lbest = l;
    }
  }

  // translate^{lbest}(s) == best (== rep), so s == translate^{N-lbest}(rep)
  *l_out = (N - lbest) % N;

  return best;
}

static int compute_period(uint64_t rep, int N, uint64_t mask) {
  uint64_t t = rep;

  for (int r = 1; r <= N; r++) {
    t = translate(t, N, mask);
    if (t == rep) {
      return r;
    }
  }

  return N; // unreachable: translate^N is always identity
}

static inline int popcount64(uint64_t x) {
  int c = 0;
  while (x) {
    x &= x - 1;
    c++;
  }

  return c;
}

/* sector construction */

spin_sector_t *spin_sector_build(int N, int nup, int k) {
  if (N <= 0 || N > 63 || nup < 0 || nup > N || k < 0 || k >= N) {
    return NULL;
  }

  uint64_t mask = (1ULL << N) - 1ULL;
  uint64_t dim_full = 1ULL << N;

  uint64_t *reps = NULL;
  int *periods = NULL;
  int count = 0, cap = 0;

  for (uint64_t s = 0; s < dim_full; s++) {
    if (popcount64(s) != nup) {
      continue;
    }

    int l;
    uint64_t rep = find_representative(s, N, mask, &l);
    if (rep != s) {
      continue; // s is not the representative of its own orbit
    }

    int R = compute_period(rep, N, mask);
    if ((k * R) % N != 0) {
      continue; // incompatible: no nonzero-norm momentum state on this orbit
    }

    if (count == cap) {
      cap = cap ? cap * 2 : 64;

      uint64_t *nr = realloc(reps, (size_t)cap * sizeof *reps);
      int *np = realloc(periods, (size_t)cap * sizeof *periods);
      if (!nr || !np) {
        free(nr ? nr : reps);
        free(np ? np : periods);

        return NULL;
      }

      reps = nr;
      periods = np;
    }

    reps[count] = rep;
    periods[count] = R;
    count++;
  }

  spin_sector_t *sec = malloc(sizeof *sec);
  if (!sec) {
    free(reps);
    free(periods);

    return NULL;
  }

  sec->N = N;
  sec->nup = nup;
  sec->k = k;
  sec->dim = count;
  sec->rep = reps;
  sec->period = periods;

  return sec;
}

void spin_sector_free(spin_sector_t *sec) {
  if (!sec) {
    return;
  }

  free(sec->rep);
  free(sec->period);
  free(sec);
}

int spin_sector_find(const spin_sector_t *sec, uint64_t rep) {
  if (!sec) {
    return -1;
  }

  int lo = 0, hi = sec->dim - 1;
  while (lo <= hi) {
    int mid = lo + (hi - lo) / 2;
    if (sec->rep[mid] == rep) {
      return mid;
    }

    if (sec->rep[mid] < rep) {
      lo = mid + 1;
    } else {
      hi = mid - 1;
    }
  }

  return -1;
}

/* matrix-free Hamiltonian in one symmetry sector
 *
 * For representative row a, walking every bond once gives H|rep_a> as a short
 * list of (image state, coefficient) pairs - at most N+1 entries  (1 diagonal
 * Jz accumulator + up to N flip-flop terms). Each image state maps to some
 * (rep_b, l) via find_representative, and the sector matrix element is :
 *
 *   <k,b|H|k,a> += c_i * \sqrt(R_a / R_b) * \exp(-i * 2 * \pi * k * l / N)
 *
 * INFO: Derivation :
 *  expand |k,a> = (1 / \sqrt(R_a)) \sum_p \exp(i * 2 * \pi * k * p  / N)
 *   T^p|a>, commute H past each T^p onto the \bra, and use T|k,b> =
 *  \exp(-i * 2 * \pi * k / N)|k,b> to collapse resulting sum over p into a
 *  factor of R_a.
 *  Remains is \sqrt(R_a) * <k,b|H|a>, and expanding only \bra  <k,b| = (1 /
 *  \sqrt(R_b)) \sum_p' \exp(-i * 2 * \pi*k * p' / N) <T^p' b| against (already
 *  computational-basis) image states from H|a> gives \sqrt(R_a / R_b) shown
 *  above.
 *
 * Row a's nonzero list depends only on rep[a] so both nnz-counting pass and
 * fill pass are row-parallel: give each row its own pre-computed disjoint slice
 * of col_ind/values and there is nothing to lock, no scatter/gather race, and
 * no thread-private buffer to reduce afterward.
 */
sparse_matrix_t *spin_sector_hamiltonian(const spin_sector_t *sec, double Jxy,
                                         double Jz, int pbc) {
  if (!sec || sec->dim <= 0) {
    return NULL;
  }

  int N = sec->N, dim = sec->dim;
  uint64_t mask = (1ULL << N) - 1ULL;
  int nbonds = pbc ? N : N - 1;

  int *nnz_row = malloc((size_t)dim * sizeof *nnz_row);
  if (!nnz_row) {
    return NULL;
  }

  /* Pass 1: count nonzeros per row (1 diagonal slot + one slot per
   * spin-flipping bond). Purely a function of rep[a]; independent rows. */
#pragma omp parallel for schedule(guided)
  for (int a = 0; a < dim; a++) {
    uint64_t s = sec->rep[a];
    int cnt = 1; // diagonal

    for (int b = 0; b < nbonds; b++) {
      int j = b, kk = (b + 1) % N;

      if (((s >> j) & 1ULL) != ((s >> kk) & 1ULL)) {
        cnt++;
      }
    }

    nnz_row[a] = cnt;
  }

  int *row_ptr = malloc((size_t)(dim + 1) * sizeof *row_ptr);
  if (!row_ptr) {
    free(nnz_row);

    return NULL;
  }

  row_ptr[0] = 0;
  for (int a = 0; a < dim; a++) {
    row_ptr[a + 1] = row_ptr[a] + nnz_row[a];
  }

  int nnz = row_ptr[dim];

  free(nnz_row);

  sparse_matrix_t *H = sparse_alloc(dim, dim, nnz);
  if (!H) {
    free(row_ptr);

    return NULL;
  }

  for (int a = 0; a <= dim; a++) {
    H->row_ptr[a] = row_ptr[a];
  }

  free(row_ptr);

  /* Pass 2: fill each row's already-known, disjoint slice in parallel. */
#pragma omp parallel for schedule(guided)
  for (int a = 0; a < dim; a++) {
    uint64_t s = sec->rep[a];
    int Ra = sec->period[a];
    int idx = H->row_ptr[a];

    double diag = 0.0;
    for (int b = 0; b < nbonds; b++) {
      int j = b;
      int kk = (b + 1) % N;
      diag += Jz * sz_bit(s, j) * sz_bit(s, kk);
    }

    H->col_ind[idx] = a;
    H->values[idx] = c_new(diag, 0.0);
    idx++;

    for (int b = 0; b < nbonds; b++) {
      int j = b, kk = (b + 1) % N;
      int bj = (int)((s >> j) & 1ULL), bk = (int)((s >> kk) & 1ULL);
      if (bj == bk) {
        continue;
      }

      uint64_t flipped = s ^ (1ULL << j) ^ (1ULL << kk);
      int l;
      uint64_t rep_b = find_representative(flipped, N, mask, &l);
      int bcol = spin_sector_find(sec, rep_b);
      if (bcol < 0) {
        /* NOTE: H commutes with T, so any state reachable from a compatible
         * representative is itself compatible. Emit an explicit zero rather
         * than leaving the slot uninitialized. */
        H->col_ind[idx] = a;
        H->values[idx] = c_zero();
        idx++;

        continue;
      }

      int Rb = sec->period[bcol];
      double amp = 0.5 * Jxy * sqrt((double)Ra / (double)Rb);
      double phase = -2.0 * M_PI * (double)sec->k * (double)l / (double)N;

      H->col_ind[idx] = bcol;
      H->values[idx] = c_scale(c_new(cos(phase), sin(phase)), amp);
      idx++;
    }
  }

  return H;
}

/* S^z_q excitation */

int spin_apply_szq(const spin_sector_t *src, const cvector_t *psi0, int q_index,
                   spin_sector_t **out_target, cvector_t **out_phi0,
                   double *out_I0) {
  if (!src || !psi0 || psi0->n != src->dim || !out_target || !out_phi0 ||
      !out_I0) {
    return 1;
  }

  int N = src->N;
  int q_mod = ((q_index % N) + N) % N;
  int k_target = (src->k + q_mod) % N;

  spin_sector_t *target = spin_sector_build(N, src->nup, k_target);
  if (!target) {
    return 1;
  }

  cvector_t *phi0 = cvector_alloc(target->dim);
  if (!phi0) {
    spin_sector_free(target);

    return 1;
  }
  for (int i = 0; i < target->dim; i++) {
    phi0->data[i] = c_zero();
  }

  double q = 2.0 * M_PI * (double)q_mod / (double)N;

  for (int a = 0; a < src->dim; a++) {
    complex_t ca = psi0->data[a];
    if (ca.re == 0.0 && ca.im == 0.0) {
      continue;
    }

    uint64_t s = src->rep[a];
    /* NOTE: c_a = \sum_j \exp(i * q * j) * S^z_j(s); this vanishes
     * automatically when q_index*R_a is not a multiple of N, which is the
     * compatibility condition for rep a to also carry a valid k_target state.
     */
    complex_t weight = c_zero();
    for (int j = 0; j < N; j++) {
      double phase = q * (double)j;

      complex_t ej = c_new(cos(phase), sin(phase));
      weight = c_add(weight, c_scale(ej, sz_bit(s, j)));
    }

    if (c_abs(weight) < 1e-12) {
      continue;
    }

    int bcol = spin_sector_find(target, s);
    if (bcol < 0) {
      continue; // weight should be ~0 here; skip defensively
    }

    complex_t contrib = c_mul(ca, weight);
    phi0->data[bcol] = c_add(phi0->data[bcol], contrib);
  }

  double I0 = 0.0;
  for (int i = 0; i < target->dim; i++) {
    I0 += c_abs2(phi0->data[i]);
  }

  *out_target = target;
  *out_phi0 = phi0;
  *out_I0 = I0;

  return 0;
}

/* continued-fraction DSF */

double spin_dsf_continued_fraction(const double *alpha, const double *beta,
                                   int m, double E0, double I0, double omega,
                                   double eta) {
  if (!alpha || m < 1 || I0 < 0.0 || eta <= 0.0) {
    return 0.0;
  }

  complex_t z = c_new(E0 + omega, eta);
  complex_t one = c_new(1.0, 0.0);

  /* Bottom-up evaluation: g_{m-1} = 1/(z - alpha_{m-1}), then
   * g_i = 1/(z - alpha_i - beta_i^2 * g_{i+1}) for i = m-2 downto 0. */
  complex_t g = c_div(one, c_sub(z, c_new(alpha[m - 1], 0.0)));

  for (int i = m - 2; i >= 0; i--) {
    complex_t denom =
        c_sub(c_sub(z, c_new(alpha[i], 0.0)), c_scale(g, beta[i] * beta[i]));
    g = c_div(one, denom);
  }

  complex_t G = c_scale(g, I0);

  return -G.im / M_PI;
}

/* Reflection (spatial parity) symmetry */

/* Reflection R: site j -> N-1-j, applied directly to a computational basis
 * bitstring.
 * NOTE: this is a different operation from spin inversion (bit-flip s -> s ^
 * mask, i.e. S^z_j -> -S^z_j on every site) - two symmetries are independent
 * and must not be conflated.
 */
static uint64_t reflect_bits(uint64_t s, int N) {
  uint64_t r = 0;
  for (int j = 0; j < N; j++) {
    if ((s >> j) & 1ULL) {
      r |= 1ULL << (N - 1 - j);
    }
  }
  return r;
}

cmatrix_t *spin_sector_reflection_matrix(const spin_sector_t *sec) {
  if (!sec) {
    return NULL;
  }

  int N = sec->N, dim = sec->dim, k = sec->k;
  int is_zero = (k == 0);
  int is_pi = (N % 2 == 0) && (k == N / 2);
  if (!is_zero && !is_pi) {
    return NULL; /* no reflection quantum number at this momentum */
  }

  uint64_t mask = (N == 64) ? ~0ULL : ((1ULL << N) - 1ULL);
  double theta = 2.0 * M_PI * (double)k / (double)N;

  cmatrix_t *R = cmatrix_alloc(dim, dim);
  if (!R) {
    return NULL;
  }
  for (int i = 0; i < dim * dim; i++) {
    R->data[i] = c_zero();
  }

  /* Derivation: writing normalized momentum state as :
   *  |k,a> = (\sqrt(R_a) / N) * \sum_{r=0}^{N-1} \exp(-i * \theta * r) T^r
   *            |rep_a>
   * and using dihedral relation R T^r = T^{-r} R, one gets :
   *  R|k,a> = (\sqrt(R_a) / N) * \exp(-i * \theta * l0) * \sum_r \exp(+i *
   *             \theta * r) T^r |rep_b>
   *
   * Where
   *   rep_b = representative of s0 = R|rep_a>
   *   T^{l0}(rep_b) = s0
   *  At k=0 or k=N/2, \theta * r is always an integer multiple of \pi, so
   * \exp(+i * \theta * r) == \exp(-i* \theta * r) exactly - this is why a
   * reflection quantum number only exists at these two momenta and sum
   * collapses to (N / \sqrt(R_b)) |k,b>, giving <k,b|R|k,a> = \sqrt(R_a / R_b)
   * * \exp(-i * \theta * l0).
   */
  for (int a = 0; a < dim; a++) {
    uint64_t s0 = reflect_bits(sec->rep[a], N);
    int l0;
    uint64_t rep_b = find_representative(s0, N, mask, &l0);
    int b = spin_sector_find(sec, rep_b);
    if (b < 0) {
      /* Should not happen: H (and hence sector's compatibility condition) is
       * reflection-symmetric, so a compatible rep's mirror image is always
       * itself compatible at k=0/N/2. */
      cmatrix_free(R);

      return NULL;
    }

    double Ra = (double)sec->period[a];
    double Rb = (double)sec->period[b];
    double phase = -theta * (double)l0;
    complex_t val = c_scale(c_new(cos(phase), sin(phase)), sqrt(Ra / Rb));
    CMAT(R, b, a) = val;
  }

  return R;
}

cmatrix_t *spin_sector_parity_project(const spin_sector_t *sec,
                                      const sparse_matrix_t *H, int parity,
                                      int *dim_out) {
  if (!sec || !H || (parity != 1 && parity != -1) || !dim_out) {
    return NULL;
  }

  cmatrix_t *R = spin_sector_reflection_matrix(sec);
  if (!R) {
    return NULL;
  }

  int dim = sec->dim;

  cmatrix_t *P = cmatrix_alloc(dim, dim);
  if (!P) {
    cmatrix_free(R);

    return NULL;
  }

  double half_parity = 0.5 * (double)parity;
  double trace = 0.0;

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      complex_t rij = CMAT(R, i, j);
      complex_t pij = c_scale(rij, half_parity);

      if (i == j) {
        pij = c_add(pij, c_new(0.5, 0.0));
      }

      CMAT(P, i, j) = pij;
    }

    trace += CMAT(P, i, i).re;
  }

  cmatrix_free(R);

  int expected_dim = (int)(trace + 0.5); // Tr(P) = rank(P) for a projector
  *dim_out = expected_dim;

  if (expected_dim <= 0) {
    cmatrix_free(P);

    return cmatrix_alloc(0, 0);
  }

  /* NOTE: Modified Gram-Schmidt with pivoting on P's own columns, stopping once
   * expected_dim orthonormal vectors are found. */
  cmatrix_t *V = cmatrix_alloc(dim, expected_dim);
  if (!V) {
    cmatrix_free(P);

    return NULL;
  }

  complex_t *col = malloc((size_t)dim * sizeof *col);
  int accepted = 0;
  const double norm_tol = 1e-8;
  for (int c = 0; c < dim && accepted < expected_dim; c++) {
    for (int i = 0; i < dim; i++) {
      col[i] = CMAT(P, i, c);
    }

    for (int a = 0; a < accepted; a++) {
      complex_t dot = c_zero();

      for (int i = 0; i < dim; i++) {
        dot = c_add(dot, c_mul(c_conj(CMAT(V, i, a)), col[i]));
      }

      for (int i = 0; i < dim; i++) {
        col[i] = c_sub(col[i], c_mul(dot, CMAT(V, i, a)));
      }
    }

    double norm2 = 0.0;
    for (int i = 0; i < dim; i++) {
      norm2 += c_abs2(col[i]);
    }

    double norm = sqrt(norm2);
    if (norm < norm_tol) {
      continue; // this column was (numerically) already in span
    }

    for (int i = 0; i < dim; i++) {
      CMAT(V, i, accepted) = c_scale(col[i], 1.0 / norm);
    }

    accepted++;
  }

  free(col);
  cmatrix_free(P);

  if (accepted != expected_dim) {
    cmatrix_free(V);

    return NULL;
  }

  cmatrix_t *Hd = cmatrix_alloc(dim, dim);
  if (!Hd) {
    cmatrix_free(V);

    return NULL;
  }

  for (int i = 0; i < dim * dim; i++) {
    Hd->data[i] = c_zero();
  }

  for (int i = 0; i < dim; i++) {
    for (int idx = H->row_ptr[i]; idx < H->row_ptr[i + 1]; idx++) {
      int j = H->col_ind[idx];

      CMAT(Hd, i, j) = c_add(CMAT(Hd, i, j), H->values[idx]);
    }
  }

  cmatrix_t *Vd = cmatrix_adjoint(V);
  cmatrix_t *HV = cmatrix_multiply(Hd, V);
  cmatrix_t *Hblock = (Vd && HV) ? cmatrix_multiply(Vd, HV) : NULL;

  cmatrix_free(V);
  cmatrix_free(Vd);
  cmatrix_free(HV);
  cmatrix_free(Hd);

  return Hblock;
}
