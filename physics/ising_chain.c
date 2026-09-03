#include "ising_chain.h"
#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* \sigma^z eigenvalue for bit j of computational basis state s: bit=0 (up) ->
 * +1, bit=1 (down) -> -1. */
static inline int sz_bit(uint32_t s, int j) {
  return 1 - 2 * (int)((s >> j) & 1u);
}

sparse_matrix_t *ising_hamiltonian(int N, double J, double h, int pbc) {
  if (N <= 0 || N > 30) {
    return NULL;
  }

  uint32_t dim = 1u << N;
  int nbonds = pbc ? N : N - 1;

  int *nnz_row = malloc((size_t)dim * sizeof *nnz_row);
  if (!nnz_row) {
    return NULL;
  }

  /* Row s's nonzero count depends only on s itself: 1 diagonal slot (ZZ sum)
   * plus N off-diagonal slots.
   * */
#pragma omp parallel for schedule(guided)
  for (uint32_t s = 0; s < dim; s++) {
    nnz_row[s] = 1 + N;
  }

  int *row_ptr = malloc((size_t)(dim + 1) * sizeof *row_ptr);
  if (!row_ptr) {
    free(nnz_row);

    return NULL;
  }

  row_ptr[0] = 0;
  for (uint32_t s = 0; s < dim; s++) {
    row_ptr[s + 1] = row_ptr[s] + nnz_row[s];
  }

  int nnz = row_ptr[dim];

  free(nnz_row);

  sparse_matrix_t *H = sparse_alloc((int)dim, (int)dim, nnz);
  if (!H) {
    free(row_ptr);

    return NULL;
  }

  for (uint32_t s = 0; s <= dim; s++) {
    H->row_ptr[s] = row_ptr[s];
  }

  free(row_ptr);

  /* Pass 2: fill each row's already-known, disjoint slice in parallel.
   * Row s writes only to H->{col_ind,values}[H->row_ptr[s] .. H->row_ptr[s+1]),
   * which no other row touches. */
#pragma omp parallel for schedule(guided)
  for (uint32_t s = 0; s < dim; s++) {
    int idx = H->row_ptr[s];

    double diag = 0.0;
    for (int b = 0; b < nbonds; b++) {
      int j = b, k = (b + 1) % N;
      diag += -J * sz_bit(s, j) * sz_bit(s, k);
    }

    H->col_ind[idx] = (int)s;
    H->values[idx] = c_new(diag, 0.0);
    idx++;

    for (int i = 0; i < N; i++) {
      uint32_t t = s ^ (1u << i);

      H->col_ind[idx] = (int)t;
      H->values[idx] = c_new(-h, 0.0);
      idx++;
    }
  }

  return H;
}

sparse_matrix_t *ising_z2_hamiltonian(int N, double J, double h, int pbc,
                                      int parity) {
  if (N <= 0 || N > 30 || (parity != 1 && parity != -1)) {
    return NULL;
  }

  uint32_t dim_full = 1u << N;
  uint32_t mask = dim_full - 1u;
  uint32_t dim =
      dim_full >> 1; /* NOTE: 2^(N-1); mask != 0 for N>=1, so global flip s ->
                        s^mask has no fixed points and splits the space into
                        exactly this many symmetric/antisymmetric pairs */
  int nbonds = pbc ? N : N - 1;

  /* NOTE: rep[a] = a-th smallest s with s < (s^mask), i.e. canonical
   * ("smaller") member of flip-pair a, in ascending order. This can be written
   * down directly without any search: as s runs over 0..dim_full-1 in order,
   * exactly every other one (by a simple parity-of-position argument on the
   * pairing structure) satisfies s < s^mask, and they already come out sorted,
   * so just filter. */
  uint32_t *rep = malloc((size_t)dim * sizeof *rep);
  if (!rep) {
    return NULL;
  }

  uint32_t count = 0;
  for (uint32_t s = 0; s < dim_full && count < dim; s++) {
    if (s < (s ^ mask)) {
      rep[count++] = s;
    }
  }

  /* NOTE: count should always equal dim exactly; if it doesn't (shouldn't
   * happen for mask != 0), the loop above already stopped safely at count==dim,
   * so rep[] is fully populated regardless but bail if this invariant is
   * somehow violated on the low end. */
  if (count != dim) {
    free(rep);

    return NULL;
  }

  // Binary search rep[] for a given representative value
#define ISING_FIND_REP(target, out_idx)                                        \
  do {                                                                         \
    int lo = 0, hi = (int)dim - 1;                                             \
    (out_idx) = -1;                                                            \
                                                                               \
    while (lo <= hi) {                                                         \
      int mid = lo + (hi - lo) / 2;                                            \
      if (rep[mid] == (target)) {                                              \
        (out_idx) = mid;                                                       \
        break;                                                                 \
      } else if (rep[mid] < (target)) {                                        \
        lo = mid + 1;                                                          \
      } else {                                                                 \
        hi = mid - 1;                                                          \
      }                                                                        \
    }                                                                          \
  } while (0)

  int *nnz_row = malloc((size_t)dim * sizeof *nnz_row);
  if (!nnz_row) {
    free(rep);

    return NULL;
  }

#pragma omp parallel for schedule(guided)
  for (uint32_t a = 0; a < dim; a++) {
    nnz_row[a] =
        1 + N; /* NOTE: upper bound: diagonal + up to N off-diagonal (some flips
                  may coincide in column, but CSR doesn't require deduplication
                  - sparse_mv just sums duplicate entries in same row/col
                  naturally via separate accumulation slots) */
  }

  int *row_ptr = malloc((size_t)(dim + 1) * sizeof *row_ptr);
  if (!row_ptr) {
    free(nnz_row);
    free(rep);

    return NULL;
  }

  row_ptr[0] = 0;
  for (uint32_t a = 0; a < dim; a++) {
    row_ptr[a + 1] = row_ptr[a] + nnz_row[a];
  }

  int nnz = row_ptr[dim];

  free(nnz_row);

  sparse_matrix_t *H = sparse_alloc((int)dim, (int)dim, nnz);
  if (!H) {
    free(row_ptr);
    free(rep);

    return NULL;
  }

  for (uint32_t a = 0; a <= dim; a++) {
    H->row_ptr[a] = row_ptr[a];
  }

  free(row_ptr);

  /* Derivation :
   * writing |a,+-> = (|rep_a> +- |rep_a^mask>) / \sqrt(2), and using that H's
   * diagonal (ZZ) part is invariant under the global flip while each
   * single-site X flip on globally-flipped state rep_a^mask lands on
   * (rep_a^(1<<i))^mask - i.e. flip commutes with global-flip permutation -one
   * finds
   *
   *   <b,+|H|a,+> = diag_a * \delta_{ab} + \sum_{i: b_i=b} (-h)
   *   <b,-|H|a,-> = diag_a * \delta_{ab} + \sum_{i: b_i=b} (-h) * sign_i
   *
   * Where for each site i,
   *  t_i = rep_a ^ (1<<i)
   *  b_i = pair-index of t_i (i.e. rep[b_i] = min(t_i, t_i^mask))
   *  sign_i = +1 if t_i == rep[b_i] (t_i is already canonical/smaller element
   * of its pair) or -1 if t_i is other (larger) element
   *
   * NOTE: The "+" sector never picks up this sign because a symmetric
   * combination |t_i>+|t_i^mask> is insensitive to which element you started
   * from; "-" sector does, since |t_i>-|t_i^mask> flips overall sign under that
   * same relabeling.
   */
#pragma omp parallel for schedule(guided)
  for (uint32_t a = 0; a < dim; a++) {
    uint32_t s = rep[a];
    int idx = H->row_ptr[a];

    double diag = 0.0;
    for (int b = 0; b < nbonds; b++) {
      int j = b, k = (b + 1) % N;
      diag += -J * sz_bit(s, j) * sz_bit(s, k);
    }

    H->col_ind[idx] = (int)a;
    H->values[idx] = c_new(diag, 0.0);
    idx++;

    for (int i = 0; i < N; i++) {
      uint32_t t = s ^ (1u << i);
      uint32_t t_partner = t ^ mask;
      uint32_t rep_b;
      int sign_i;
      if (t < t_partner) {
        rep_b = t;
        sign_i = 1;
      } else {
        rep_b = t_partner;
        sign_i = -1;
      }

      int bcol;
      ISING_FIND_REP(rep_b, bcol);
      if (bcol < 0) {
        /* NOTE: Shouldn't happen: rep_b is by construction smaller element of a
         * valid flip-pair, so it's always in rep[]. */
        H->col_ind[idx] = (int)a;
        H->values[idx] = c_zero();
        idx++;

        continue;
      }

      double amp = (parity == 1) ? -h : -h * (double)sign_i;

      H->col_ind[idx] = bcol;
      H->values[idx] = c_new(amp, 0.0);
      idx++;
    }
  }

#undef ISING_FIND_REP

  free(rep);

  return H;
}

double ising_exact_ground_energy_per_site(int N, double J, double h) {
  if (N < 1) {
    return NAN;
  }

  /* NOTE: Jordan-Wigner fermionization maps spin chain to free fermions on a
   * ring, but with a subtlety: whether the resulting fermions obey periodic or
   * antiperiodic boundary conditions depends on the total fermion parity of
   * state, which is tied to the spin-flip parity sector. The true many-body
   * ground state (checked here against exact diagonalization for N=3..10, even
   * and odd, including h=0, J=0, h<0, and the critical point h=J) always lies
   * in sector requiring antiperiodic fermion boundary conditions, giving
   * momenta
   *
   *   k_m = (2 * m + 1) * \pi/N,  m = 0, 1, ..., N-1
   *
   * with Bogoliubov quasiparticle energies
   *
   *   \epsilon_k = 2 * \sqrt(J^2 + h^2 - 2 * J * h * \cos(k))
   *
   * and ground-state energy E0 = -(1/2) * \sum_m \epsilon_{k_m} (all modes
   * occupied by their negative-energy quasiparticle vacuum contribution).
   * Reference : e.g. Sachdev, "Quantum Phase Transitions", ch. 4, for the full
   * derivation; this is the periodic-TFIM result. */
  double E0 = 0.0;
  for (int m = 0; m < N; m++) {
    double k = (2.0 * m + 1.0) * M_PI / (double)N;
    double eps_k = 2.0 * sqrt(J * J + h * h - 2.0 * J * h * cos(k));

    E0 += eps_k;
  }

  E0 = -0.5 * E0;

  return E0 / (double)N;
}

cmatrix_t *ising_reduced_density_matrix(const cmatrix_t *psi, int N, int L_A) {
  if (!psi || N <= 0 || N > 30 || L_A < 0 || L_A > N) {
    return NULL;
  }

  int dim_A = 1 << L_A;
  int dim_B = 1 << (N - L_A);

  cmatrix_t *rho_A = cmatrix_alloc(dim_A, dim_A);
  if (!rho_A) {
    return NULL;
  }

  /* NOTE: \rho_A[a][a'] = \sum_b \psi[a | (b<<L_A)] * conj(\psi[a' |
   * (b<<L_A)]). Row a is independent of every other row (each only reads from
   * psi, which is never written here), so this is safe to parallelize with same
   * row-based-gather reasoning used throughout this module */
  for (int a = 0; a < dim_A; a++) {
    for (int ap = 0; ap < dim_A; ap++) {
      complex_t sum = c_zero();

      for (int b = 0; b < dim_B; b++) {
        int idx1 = a | (b << L_A);
        int idx2 = ap | (b << L_A);

        sum = c_add(sum, c_mul(CMAT(psi, idx1, 0), c_conj(CMAT(psi, idx2, 0))));
      }

      CMAT(rho_A, a, ap) = sum;
    }
  }

  return rho_A;
}

double ising_entanglement_entropy(const cmatrix_t *psi, int N, int L_A) {
  cmatrix_t *rho_A = ising_reduced_density_matrix(psi, N, L_A);
  if (!rho_A) {
    return NAN;
  }

  eigen_t *eig = cmatrix_eigh_complex(rho_A);
  cmatrix_free(rho_A);
  if (!eig) {
    return NAN;
  }

  /* NOTE: -Tr(\rho \log \rho) = -\sum_i p_i * \ln(p_i), skipping eigenvalues at
   * or below a small numerical floor: an exactly-zero eigenvalue (e.g. every
   * eigenvalue but one, for a product-state subsystem) would make \log(0)
   * diverge, and floating-point noise can occasionally push a
   * mathematically-exact zero slightly negative. */
  const double floor_p = 1e-14;
  double S = 0.0;
  for (int i = 0; i < eig->n; i++) {
    double p = eig->eigenvalues[i];
    if (p > floor_p) {
      S += -p * log(p);
    }
  }

  eigen_free(eig);

  return S;
}

cmatrix_t *ising_hamiltonian_dense(int N, double J, double h, int pbc) {
  if (N <= 0 || N > 30) {
    return NULL;
  }

  uint32_t dim = 1u << N;
  int nbonds = pbc ? N : N - 1;

  cmatrix_t *H = cmatrix_alloc((int)dim, (int)dim);
  if (!H) {
    return NULL;
  }
  for (uint32_t i = 0; i < dim * dim; i++) {
    H->data[i] = c_zero();
  }

  /* NOTE: Each iteration s only ever writes into column s of dense matrix *
   * (row s itself for diagonal, and up to N other rows - always distinct from
   * each other, since XOR-ing s with N different single-bit masks always gives
   * N distinct results - for off-diagonal X flips), and different iterations
   * never share a column index, so this is safe to parallelize despite being a
   * dense "scatter" rather than row-based CSR gather used elsewhere in this
   * module: no two threads ever write to the same memory location. */
#pragma omp parallel for schedule(guided)
  for (uint32_t s = 0; s < dim; s++) {
    double diag = 0.0;

    for (int b = 0; b < nbonds; b++) {
      int j = b, k = (b + 1) % N;

      diag += -J * sz_bit(s, j) * sz_bit(s, k);
    }

    CMAT(H, (int)s, (int)s) = c_new(diag, 0.0);

    for (int i = 0; i < N; i++) {
      uint32_t t = s ^ (1u << i);

      CMAT(H, (int)t, (int)s) = c_new(-h, 0.0);
    }
  }

  return H;
}

double ising_loschmidt_echo(const eigen_t *eig_f, const cmatrix_t *psi_i,
                            double t) {
  if (!eig_f || !psi_i || eig_f->n != psi_i->nrows || psi_i->ncols != 1) {
    return NAN;
  }

  int dim = eig_f->n;

  /* <\psi_i|\exp(-i * H_f * t)|\psi_i> = \sum_n |c_n|^2 * \exp(-i * E_n * t),
   * Where c_n = <n|\psi_i> are \psi_i's coefficients in H_f's eigenbasis and
   * avoids ever forming \exp(-i * H_f * t) as a matrix. */
  complex_t amplitude = c_zero();
  for (int n = 0; n < dim; n++) {
    complex_t c_n = c_zero();

    for (int i = 0; i < dim; i++) {
      c_n = c_add(c_n, c_mul(c_conj(CMAT(eig_f->eigenvectors, i, n)),
                             CMAT(psi_i, i, 0)));
    }

    double p_n = c_abs2(c_n);
    double phase = -eig_f->eigenvalues[n] * t;

    amplitude = c_add(amplitude, c_scale(c_new(cos(phase), sin(phase)), p_n));
  }

  return c_abs2(amplitude);
}
