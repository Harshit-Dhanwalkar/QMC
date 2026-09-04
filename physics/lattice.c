/*
Tight-binding lattice models: dense Hamiltonian builders + exact analytic
dispersions (Bloch's theorem / standing-wave quantization) for validation.
*/

#include "lattice.h"
#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "../core/matrix.h"
#include "../core/random.h"
#include "../core/vector.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

cmatrix_t *lattice_build_1d_chain(int n_sites, double epsilon0, double t,
                                  lattice_bc_t bc) {
  if (n_sites < 1) {
    return NULL;
  }

  cmatrix_t *H = cmatrix_alloc(n_sites, n_sites);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < n_sites; i++) {
    for (int j = 0; j < n_sites; j++) {
      CMAT(H, i, j) = c_zero();
    }

    CMAT(H, i, i) = c_real(epsilon0);
  }

  for (int i = 0; i + 1 < n_sites; i++) {
    CMAT(H, i, i + 1) = c_real(-t);
    CMAT(H, i + 1, i) = c_real(-t);
  }

  if (bc == LATTICE_PERIODIC && n_sites > 1) {
    CMAT(H, 0, n_sites - 1) = c_add(CMAT(H, 0, n_sites - 1), c_real(-t));
    CMAT(H, n_sites - 1, 0) = c_add(CMAT(H, n_sites - 1, 0), c_real(-t));
  }

  return H;
}

/* \cos(angle_n) for n-th mode of an n_sites-site chain with boundary condition
 * bc; shared per-direction piece of both 1D and 2D analytic dispersions. Writes
 * n_sites values into cos_out. */
static void chain_cos_angles(int n_sites, lattice_bc_t bc, double *cos_out) {
  if (bc == LATTICE_PERIODIC) {
    for (int n = 0; n < n_sites; n++) {
      cos_out[n] = cos(2.0 * M_PI * n / n_sites);
    }
  } else {
    for (int n = 1; n <= n_sites; n++) {
      cos_out[n - 1] = cos(M_PI * n / (n_sites + 1.0));
    }
  }
}

static int cmp_double(const void *a, const void *b) {
  double da = *(const double *)a, db = *(const double *)b;

  return (da > db) - (da < db);
}

void lattice_1d_chain_analytic(int n_sites, double epsilon0, double t,
                               lattice_bc_t bc, double *E_out) {
  if (n_sites < 1 || !E_out) {
    return;
  }

  chain_cos_angles(n_sites, bc, E_out);
  for (int n = 0; n < n_sites; n++) {
    E_out[n] = epsilon0 - 2.0 * t * E_out[n];
  }

  qsort(E_out, (size_t)n_sites, sizeof(double), cmp_double);
}

cmatrix_t *lattice_build_2d_square(int nx, int ny, double epsilon0, double t,
                                   lattice_bc_t bc) {
  if (nx < 1 || ny < 1) {
    return NULL;
  }

  int N = nx * ny;
  cmatrix_t *H = cmatrix_alloc(N, N);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

#define IDX(ix, iy) ((ix) * ny + (iy))

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      int i = IDX(ix, iy);
      CMAT(H, i, i) = c_real(epsilon0);

      if (ix + 1 < nx) {
        int j = IDX(ix + 1, iy);
        CMAT(H, i, j) = c_real(-t);
        CMAT(H, j, i) = c_real(-t);
      } else if (bc == LATTICE_PERIODIC && nx > 1) {
        int j = IDX(0, iy);
        CMAT(H, i, j) = c_add(CMAT(H, i, j), c_real(-t));
        CMAT(H, j, i) = c_add(CMAT(H, j, i), c_real(-t));
      }

      if (iy + 1 < ny) {
        int j = IDX(ix, iy + 1);
        CMAT(H, i, j) = c_real(-t);
        CMAT(H, j, i) = c_real(-t);
      } else if (bc == LATTICE_PERIODIC && ny > 1) {
        int j = IDX(ix, 0);
        CMAT(H, i, j) = c_add(CMAT(H, i, j), c_real(-t));
        CMAT(H, j, i) = c_add(CMAT(H, j, i), c_real(-t));
      }
    }
  }

#undef IDX

  return H;
}

void lattice_2d_square_analytic(int nx, int ny, double epsilon0, double t,
                                lattice_bc_t bc, double *E_out) {
  if (nx < 1 || ny < 1 || !E_out) {
    return;
  }

  double *cos_x = malloc((size_t)nx * sizeof *cos_x);
  double *cos_y = malloc((size_t)ny * sizeof *cos_y);
  if (!cos_x || !cos_y) {
    free(cos_x);
    free(cos_y);

    return;
  }

  chain_cos_angles(nx, bc, cos_x);
  chain_cos_angles(ny, bc, cos_y);

  int k = 0;
  for (int m = 0; m < nx; m++) {
    for (int n = 0; n < ny; n++) {
      E_out[k++] = epsilon0 - 2.0 * t * cos_x[m] - 2.0 * t * cos_y[n];
    }
  }

  free(cos_x);
  free(cos_y);

  qsort(E_out, (size_t)(nx * ny), sizeof(double), cmp_double);
}

cmatrix_t *lattice_build_anderson_1d(int n_sites, double t, double disorder_W,
                                     uint64_t seed) {
  if (n_sites < 1 || disorder_W < 0.0) {
    return NULL;
  }

  cmatrix_t *H = cmatrix_alloc(n_sites, n_sites);
  if (!H) {
    return NULL;
  }

  rng_state_t rng;
  rng_seed(&rng, seed);

  for (int i = 0; i < n_sites; i++) {
    for (int j = 0; j < n_sites; j++) {
      CMAT(H, i, j) = c_zero();
    }
    double eps_i = rng_uniform_range(&rng, -disorder_W / 2.0, disorder_W / 2.0);
    CMAT(H, i, i) = c_real(eps_i);
  }

  for (int i = 0; i + 1 < n_sites; i++) {
    CMAT(H, i, i + 1) = c_real(-t);
    CMAT(H, i + 1, i) = c_real(-t);
  }

  return H;
}

double lattice_ipr(const cvector_t *psi) {
  if (!psi || psi->n < 1) {
    return 0.0;
  }

  double norm_sq = 0.0;
  for (int i = 0; i < psi->n; i++) {
    norm_sq += c_abs2(psi->data[i]);
  }
  if (norm_sq < 1e-300) {
    return 0.0;
  }

  double sum4 = 0.0;
  for (int i = 0; i < psi->n; i++) {
    double p = c_abs2(psi->data[i]) / norm_sq;
    sum4 += p * p;
  }

  return sum4;
}

cmatrix_t *lattice_build_ssh(int n_cells, double t1, double t2,
                             lattice_bc_t bc) {
  if (n_cells < 1) {
    return NULL;
  }

  int N = 2 * n_cells;
  cmatrix_t *H = cmatrix_alloc(N, N);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

  for (int c = 0; c < n_cells; c++) {
    int a = 2 * c, b = 2 * c + 1;
    CMAT(H, a, b) = c_real(-t1);
    CMAT(H, b, a) = c_real(-t1);

    if (c + 1 < n_cells) {
      int b_next = 2 * (c + 1);
      CMAT(H, b, b_next) = c_real(-t2);
      CMAT(H, b_next, b) = c_real(-t2);
    }
  }

  if (bc == LATTICE_PERIODIC && n_cells > 1) {
    CMAT(H, N - 1, 0) = c_add(CMAT(H, N - 1, 0), c_real(-t2));
    CMAT(H, 0, N - 1) = c_add(CMAT(H, 0, N - 1), c_real(-t2));
  }

  return H;
}

cmatrix_t *lattice_build_2d_square_magnetic(int nx, int ny, double epsilon0,
                                            double t, double alpha,
                                            lattice_bc_t bc_y) {
  if (nx < 1 || ny < 1) {
    return NULL;
  }

  int N = nx * ny;
  cmatrix_t *H = cmatrix_alloc(N, N);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

#define IDXM(ix, iy) ((ix) * ny + (iy))

  for (int ix = 0; ix < nx; ix++) {
    // Peierls phase for a y-bond at current x column :
    //   \exp(i * 2 * \pi * \alpha * ix)
    complex_t phase = c_exp(c_new(0.0, 2.0 * M_PI * alpha * (double)ix));
    complex_t hop_y = c_scale(phase, -t); // forward bond amplitude
    complex_t hop_y_conj = c_conj(hop_y); // reverse (h.c.) bond

    for (int iy = 0; iy < ny; iy++) {
      int i = IDXM(ix, iy);

      CMAT(H, i, i) = c_real(epsilon0);

      // x-bonds: unaffected by current gauge choice (real, always open in x)
      if (ix + 1 < nx) {
        int j = IDXM(ix + 1, iy);

        CMAT(H, i, j) = c_real(-t);
        CMAT(H, j, i) = c_real(-t);
      }

      // y-bonds: Peierls phase set by column ix
      if (iy + 1 < ny) {
        int j = IDXM(ix, iy + 1);

        CMAT(H, i, j) = hop_y;
        CMAT(H, j, i) = hop_y_conj;
      } else if (bc_y == LATTICE_PERIODIC && ny > 1) {
        int j = IDXM(ix, 0);

        CMAT(H, i, j) = c_add(CMAT(H, i, j), hop_y);
        CMAT(H, j, i) = c_add(CMAT(H, j, i), hop_y_conj);
      }
    }
  }

#undef IDXM

  return H;
}

double lattice_landau_level_energy(int n, double epsilon0, double t,
                                   double alpha) {
  double omega_c = 4.0 * M_PI * t * alpha;

  return epsilon0 - 4.0 * t + omega_c * ((double)n + 0.5);
}

cmatrix_t *lattice_hofstadter_bloch(double kx, double ky, int p, int q,
                                    double t) {
  if (p < 1 || q < 2 || p >= q) {
    return NULL;
  }

  cmatrix_t *H = cmatrix_alloc(q, q);
  if (!H) {
    return NULL;
  }

  for (int i = 0; i < q; i++) {
    for (int j = 0; j < q; j++) {
      CMAT(H, i, j) = c_zero();
    }
  }

  double alpha = (double)p / (double)q;

  /* NOTE: Landau-gauge magnetic unit cell: q sites stacked along y (index m =
   * 0..q-1). Diagonal on-site energy carries the "in-cell" y-dispersion with a
   * site-dependent Peierls phase 2*pi*alpha*m (Harper equation form);
   * off-diagonal (m, m+1) hopping is uniform bond amplitude -t, real for every
   * bond except the one that wraps magnetic unit cell back on itself (m=q-1 to
   * m=0), which alone carries Bloch phase \exp^{-i * kx} that stitches unit
   * cells together. */
  for (int m = 0; m < q; m++) {
    CMAT(H, m, m) = c_real(2.0 * t * cos(ky + 2.0 * M_PI * alpha * (double)m));
  }

  for (int m = 0; m < q - 1; m++) {
    complex_t hop = c_real(-t);
    CMAT(H, m, m + 1) = c_add(CMAT(H, m, m + 1), hop);
    CMAT(H, m + 1, m) = c_add(CMAT(H, m + 1, m), c_conj(hop));
  }

  complex_t wrap = c_scale(c_new(cos(kx), -sin(kx)), -t);
  CMAT(H, q - 1, 0) = c_add(CMAT(H, q - 1, 0), wrap);
  CMAT(H, 0, q - 1) = c_add(CMAT(H, 0, q - 1), c_conj(wrap));

  return H;
}

int lattice_hofstadter_chern_numbers(int p, int q, double t, int n_k,
                                     double *chern_out) {
  if (p < 1 || q < 2 || p >= q || n_k < 2 || !chern_out) {
    return 0;
  }

  // Eigenvectors on the n_k x n_k BZ grid, kx,ky in [0, 2* \pi)
  cmatrix_t ***U = malloc((size_t)n_k * sizeof(cmatrix_t **));
  if (!U) {
    return 0;
  }

  for (int i = 0; i < n_k; i++) {
    U[i] = malloc((size_t)n_k * sizeof(cmatrix_t *));

    if (!U[i]) {
      for (int ii = 0; ii < i; ii++) {
        free(U[ii]);
      }

      free(U);

      return 0;
    }
  }

  for (int i = 0; i < n_k; i++) {
    double kx = 2.0 * M_PI * (double)i / (double)n_k;

    for (int j = 0; j < n_k; j++) {
      double ky = 2.0 * M_PI * (double)j / (double)n_k;
      cmatrix_t *H = lattice_hofstadter_bloch(kx, ky, p, q, t);
      eigen_t *eig = cmatrix_eigh_complex(H);

      cmatrix_free(H);
      U[i][j] = cmatrix_copy(
          eig->eigenvectors); // columns = eigenvectors, already ascending order
      eigen_free(eig);
    }
  }

  for (int n = 0; n < q; n++) {
    chern_out[n] = 0.0;
  }

  for (int i = 0; i < n_k; i++) {
    int i1 = (i + 1) % n_k;

    for (int j = 0; j < n_k; j++) {
      int j1 = (j + 1) % n_k;

      for (int n = 0; n < q; n++) {
        /* NOTE: U1..U4: gauge-invariant link variables around one plaquette
         * (Refrence: Fukui-Hatsugai-Suzuki 2005). Each is the overlap of band-n
         * eigenvectors at neighboring grid points, normalized to unit modulus
         * (normalization makes gauge-invariant: any * per-k-point phase choice
         * in eigensolver's eigenvectors cancels out of |z|). */
        complex_t U1 = c_zero(), U2 = c_zero(), U3 = c_zero(), U4 = c_zero();
        for (int a = 0; a < q; a++) {
          U1 = c_add(U1,
                     c_mul(c_conj(CMAT(U[i][j], a, n)), CMAT(U[i1][j], a, n)));
          U2 = c_add(
              U2, c_mul(c_conj(CMAT(U[i1][j], a, n)), CMAT(U[i1][j1], a, n)));
          U3 = c_add(
              U3, c_mul(c_conj(CMAT(U[i1][j1], a, n)), CMAT(U[i][j1], a, n)));
          U4 = c_add(U4,
                     c_mul(c_conj(CMAT(U[i][j1], a, n)), CMAT(U[i][j], a, n)));
        }

        double n1 = c_abs(U1), n2 = c_abs(U2), n3 = c_abs(U3), n4 = c_abs(U4);
        if (n1 < 1e-14 || n2 < 1e-14 || n3 < 1e-14 || n4 < 1e-14) {
          continue; /* exact band touching at this plaquette corner: link
                       variable undefined here, skip (see docstring) */
        }

        complex_t prod = c_scale(U1, 1.0 / n1);
        prod = c_mul(prod, c_scale(U2, 1.0 / n2));
        prod = c_mul(prod, c_scale(U3, 1.0 / n3));
        prod = c_mul(prod, c_scale(U4, 1.0 / n4));

        double curvature =
            atan2(prod.im, prod.re); // principal branch, (-\pi, \pi]
        chern_out[n] += curvature;
      }
    }
  }

  for (int n = 0; n < q; n++) {
    chern_out[n] /= (2.0 * M_PI);
  }

  for (int i = 0; i < n_k; i++) {
    for (int j = 0; j < n_k; j++) {
      cmatrix_free(U[i][j]);
    }

    free(U[i]);
  }

  free(U);

  return 1;
}
