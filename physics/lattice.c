/*
Tight-binding lattice models: dense Hamiltonian builders + exact analytic
dispersions (Bloch's theorem / standing-wave quantization) for validation.
*/

#include "lattice.h"
#include "../core/complex.h"
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
