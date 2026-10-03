#include "tight_binding.h"
#include "../core/linalg/complex_eigh.h"
#include "complex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Thresholds for declaring an invariant ill-defined.
#define TB_GAP_TOL 1e-9
#define TB_LINK_TOL 1e-10
// Cap on stored eigenvector data (complex numbers) for grid routines.
#define TB_MAX_STORED_ELEMENTS 8000000UL

typedef struct {
  int i;
  int j;
  int n1;
  int n2;
  complex_t t;
} tb_hop_t;

struct tb_model {
  int dim;
  int n_orb;
  double *eps;
  tb_hop_t *hops;
  size_t n_hops;
  size_t cap;
};

// Model construction
tb_model_t *tb_model_alloc(int dim, int n_orb) {
  if (dim < 1 || dim > TB_MAX_DIM || n_orb < 1 || n_orb > TB_MAX_ORBITALS) {
    return NULL;
  }

  tb_model_t *m = calloc(1, sizeof *m);
  if (!m) {
    return NULL;
  }

  m->dim = dim;
  m->n_orb = n_orb;
  m->eps = calloc((size_t)n_orb, sizeof *m->eps);
  if (!m->eps) {
    free(m);

    return NULL;
  }

  return m;
}

void tb_model_free(tb_model_t *model) {
  if (!model) {
    return;
  }

  free(model->eps);
  free(model->hops);
  free(model);
}

int tb_model_dim(const tb_model_t *model) { return model ? model->dim : 0; }

int tb_model_orbitals(const tb_model_t *model) {
  return model ? model->n_orb : 0;
}

tb_status_t tb_set_onsite(tb_model_t *model, int orbital, double energy) {
  if (!model || orbital < 0 || orbital >= model->n_orb || !isfinite(energy)) {
    return TB_ERR_INVALID;
  }

  model->eps[orbital] = energy;

  return TB_OK;
}

tb_status_t tb_add_hopping(tb_model_t *model, int i, int j, int n1, int n2,
                           complex_t t) {
  if (!model || i < 0 || i >= model->n_orb || j < 0 || j >= model->n_orb ||
      abs(n1) > TB_MAX_HOP_RANGE || abs(n2) > TB_MAX_HOP_RANGE ||
      (model->dim == 1 && n2 != 0) || !isfinite(t.re) || !isfinite(t.im) ||
      (i == j && n1 == 0 && n2 == 0)) {
    return TB_ERR_INVALID;
  }

  if (model->n_hops == model->cap) {
    size_t new_cap = model->cap ? 2 * model->cap : 8;
    tb_hop_t *grown = realloc(model->hops, new_cap * sizeof *grown);
    if (!grown) {
      return TB_ERR_MEMORY;
    }

    model->hops = grown;
    model->cap = new_cap;
  }

  model->hops[model->n_hops++] = (tb_hop_t){i, j, n1, n2, t};

  return TB_OK;
}

// Bloch Hamiltonian and bands
static void tb_fill_hamiltonian(const tb_model_t *m, double k1, double k2,
                                cmatrix_t *H) {
  int n = m->n_orb;
  for (int i = 0; i < n; i++) {
    H->data[(size_t)i * n + i] = c_real(m->eps[i]);
  }

  for (size_t h = 0; h < m->n_hops; h++) {
    const tb_hop_t *hop = &m->hops[h];
    double phase = k1 * hop->n1 + (m->dim == 2 ? k2 * hop->n2 : 0.0);
    complex_t v = c_mul(hop->t, c_from_polar(1.0, phase));

    H->data[(size_t)hop->i * n + hop->j] =
        c_add(H->data[(size_t)hop->i * n + hop->j], v);
    H->data[(size_t)hop->j * n + hop->i] =
        c_add(H->data[(size_t)hop->j * n + hop->i], c_conj(v));
  }
}

cmatrix_t *tb_hamiltonian(const tb_model_t *model, double k1, double k2) {
  if (!model || !isfinite(k1) || !isfinite(k2)) {
    return NULL;
  }

  cmatrix_t *H = cmatrix_alloc(model->n_orb, model->n_orb);
  if (!H) {
    return NULL;
  }

  tb_fill_hamiltonian(model, k1, k2, H);

  return H;
}

// Diagonalize H(k); eigenvalues ascending, eigenvectors as columns.
static eigen_t *tb_eig(const tb_model_t *m, double k1, double k2) {
  cmatrix_t *H = tb_hamiltonian(m, k1, k2);
  if (!H) {
    return NULL;
  }

  eigen_t *eig = cmatrix_eigh_complex(H);

  cmatrix_free(H);

  return eig;
}

tb_status_t tb_bands(const tb_model_t *model, double k1, double k2,
                     double *energies) {
  if (!model || !energies || !isfinite(k1) || !isfinite(k2)) {
    return TB_ERR_INVALID;
  }

  eigen_t *eig = tb_eig(model, k1, k2);
  if (!eig) {
    return TB_ERR_MEMORY;
  }

  memcpy(energies, eig->eigenvalues, (size_t)model->n_orb * sizeof *energies);

  eigen_free(eig);

  return TB_OK;
}

tb_status_t tb_bands_along_line(const tb_model_t *model, double k1a, double k2a,
                                double k1b, double k2b, int n_points,
                                double *energies_out) {
  if (!model || !energies_out || n_points < 2 || !isfinite(k1a) ||
      !isfinite(k2a) || !isfinite(k1b) || !isfinite(k2b)) {
    return TB_ERR_INVALID;
  }

  for (int p = 0; p < n_points; p++) {
    double s = (double)p / (double)(n_points - 1);
    tb_status_t st =
        tb_bands(model, k1a + s * (k1b - k1a), k2a + s * (k2b - k2a),
                 energies_out + (size_t)p * model->n_orb);
    if (st != TB_OK) {
      return st;
    }
  }

  return TB_OK;
}

// Density of states
tb_status_t tb_dos(const tb_model_t *model, int n_k, double e_min, double e_max,
                   int n_e, double sigma, double *dos_out) {
  if (!model || !dos_out || n_k < 1 || n_k > TB_MAX_GRID || n_e < 2 ||
      n_e > (1 << 24) || !isfinite(e_min) || !isfinite(e_max) ||
      !(e_max > e_min) || !isfinite(sigma) || !(sigma > 0.0)) {
    return TB_ERR_INVALID;
  }

  double dE = (e_max - e_min) / (double)(n_e - 1);
  double window = 8.0 * sigma;
  int n_k2 = model->dim == 2 ? n_k : 1;
  double *E = malloc((size_t)model->n_orb * sizeof *E);
  if (!E) {
    return TB_ERR_MEMORY;
  }

  memset(dos_out, 0, (size_t)n_e * sizeof *dos_out);

  for (int a = 0; a < n_k; a++) {
    for (int b = 0; b < n_k2; b++) {
      double k1 = 2.0 * M_PI * (double)a / (double)n_k;
      double k2 = 2.0 * M_PI * (double)b / (double)n_k;

      tb_status_t st = tb_bands(model, k1, k2, E);
      if (st != TB_OK) {
        free(E);

        return st;
      }

      for (int band = 0; band < model->n_orb; band++) {
        double lo = ceil((E[band] - window - e_min) / dE);
        double hi = floor((E[band] + window - e_min) / dE);
        if (lo < 0.0) {
          lo = 0.0;
        }

        if (hi > (double)(n_e - 1)) {
          hi = (double)(n_e - 1);
        }

        for (int j = (int)lo; (double)j <= hi; j++) {
          double x = e_min + (double)j * dE - E[band];

          dos_out[j] += exp(-0.5 * x * x / (sigma * sigma));
        }
      }
    }
  }

  double norm = 1.0 / ((double)n_k * (double)n_k2 * sigma * sqrt(2.0 * M_PI));
  for (int j = 0; j < n_e; j++) {
    dos_out[j] *= norm;
  }

  free(E);

  return TB_OK;
}

// Berry phases
// Determinant of n x n row-major matrix A (destroyed), by LU with partial
// pivoting
static complex_t tb_det(complex_t *A, int n) {
  complex_t det = c_one();
  for (int col = 0; col < n; col++) {
    int piv = col;
    double best = c_abs2(A[(size_t)col * n + col]);

    for (int r = col + 1; r < n; r++) {
      double v = c_abs2(A[(size_t)r * n + col]);
      if (v > best) {
        best = v;
        piv = r;
      }
    }

    if (best == 0.0) {
      return c_zero();
    }

    if (piv != col) {
      for (int c = 0; c < n; c++) {
        complex_t tmp = A[(size_t)col * n + c];

        A[(size_t)col * n + c] = A[(size_t)piv * n + c];
        A[(size_t)piv * n + c] = tmp;
      }

      det = c_neg(det);
    }

    complex_t pv = A[(size_t)col * n + col];
    det = c_mul(det, pv);
    for (int r = col + 1; r < n; r++) {
      complex_t f = c_div(A[(size_t)r * n + col], pv);

      for (int c = col; c < n; c++) {
        A[(size_t)r * n + c] =
            c_sub(A[(size_t)r * n + c], c_mul(f, A[(size_t)col * n + c]));
      }
    }
  }

  return det;
}

// Unit-modulus link variable det<ua|ub> for two n x count eigenvector blocks
// (row-major, column c = band first + c)
// Sets *bad when overlap determinant is too small to define a phase. `work` has
// count * count entries
static complex_t tb_link(int n, int count, const complex_t *ua,
                         const complex_t *ub, complex_t *work, int *bad) {
  for (int c = 0; c < count; c++) {
    for (int d = 0; d < count; d++) {
      complex_t s = c_zero();

      for (int i = 0; i < n; i++) {
        s = c_add(s, c_mul(c_conj(ua[(size_t)i * count + c]),
                           ub[(size_t)i * count + d]));
      }

      work[(size_t)c * count + d] = s;
    }
  }

  complex_t det = tb_det(work, count);
  double mag = hypot(det.re, det.im);
  if (!(mag >= TB_LINK_TOL)) {
    *bad = 1;

    return c_one();
  }

  return c_scale(det, 1.0 / mag);
}

// Copy selected eigenvector columns at k into out (n x count, row-major) and
// report smaller of two gaps to neighboring bands
static tb_status_t tb_group_at(const tb_model_t *m, double k1, double k2,
                               int first, int count, complex_t *out,
                               double *gap) {
  eigen_t *eig = tb_eig(m, k1, k2);
  if (!eig) {
    return TB_ERR_MEMORY;
  }

  int n = m->n_orb;
  for (int i = 0; i < n; i++) {
    for (int c = 0; c < count; c++) {
      out[(size_t)i * count + c] =
          eig->eigenvectors->data[(size_t)i * n + first + c];
    }
  }

  double g = INFINITY;
  if (first > 0) {
    g = fmin(g, eig->eigenvalues[first] - eig->eigenvalues[first - 1]);
  }

  if (first + count < n) {
    g = fmin(g, eig->eigenvalues[first + count] -
                    eig->eigenvalues[first + count - 1]);
  }

  *gap = g;

  eigen_free(eig);

  return TB_OK;
}

static int tb_group_ok(const tb_model_t *m, int first, int count) {
  return m && first >= 0 && count >= 1 && first + count <= m->n_orb;
}

tb_status_t tb_berry_flux_grid(const tb_model_t *model, int first_band,
                               int n_bands, int n_k, double *flux_out,
                               double *min_gap) {
  if (!tb_group_ok(model, first_band, n_bands) || model->dim != 2 ||
      !flux_out || n_k < 2 || n_k > TB_MAX_GRID) {
    return TB_ERR_INVALID;
  }

  int n = model->n_orb;
  size_t blk = (size_t)n * (size_t)n_bands;
  size_t points = (size_t)n_k * (size_t)n_k;
  if (points * blk > TB_MAX_STORED_ELEMENTS) {
    return TB_ERR_INVALID;
  }

  complex_t *blocks = malloc(points * blk * sizeof *blocks);
  complex_t *work = malloc((size_t)n_bands * (size_t)n_bands * sizeof *work);
  if (!blocks || !work) {
    free(blocks);
    free(work);

    return TB_ERR_MEMORY;
  }

  double gap_min = INFINITY;
  for (int a = 0; a < n_k; a++) {
    for (int b = 0; b < n_k; b++) {
      double gap = INFINITY;
      tb_status_t st =
          tb_group_at(model, 2.0 * M_PI * (double)a / (double)n_k,
                      2.0 * M_PI * (double)b / (double)n_k, first_band, n_bands,
                      blocks + ((size_t)a * n_k + b) * blk, &gap);
      if (st != TB_OK) {
        free(blocks);
        free(work);

        return st;
      }

      gap_min = fmin(gap_min, gap);
    }
  }

  if (min_gap) {
    *min_gap = gap_min;
  }

  if (!(gap_min >= TB_GAP_TOL)) {
    free(blocks);
    free(work);

    return TB_ERR_DEGENERATE;
  }

  int bad = 0;
  for (int a = 0; a < n_k; a++) {
    int a1 = (a + 1) % n_k;
    for (int b = 0; b < n_k; b++) {
      int b1 = (b + 1) % n_k;

      const complex_t *u00 = blocks + ((size_t)a * n_k + b) * blk;
      const complex_t *u10 = blocks + ((size_t)a1 * n_k + b) * blk;
      const complex_t *u11 = blocks + ((size_t)a1 * n_k + b1) * blk;
      const complex_t *u01 = blocks + ((size_t)a * n_k + b1) * blk;
      complex_t p = tb_link(n, n_bands, u00, u10, work, &bad);

      p = c_mul(p, tb_link(n, n_bands, u10, u11, work, &bad));
      p = c_mul(p, tb_link(n, n_bands, u11, u01, work, &bad));
      p = c_mul(p, tb_link(n, n_bands, u01, u00, work, &bad));

      flux_out[(size_t)a * n_k + b] = atan2(p.im, p.re);
    }
  }

  free(blocks);
  free(work);

  return bad ? TB_ERR_DEGENERATE : TB_OK;
}

tb_status_t tb_chern_number(const tb_model_t *model, int first_band,
                            int n_bands, int n_k, double *chern_out) {
  if (!chern_out || n_k < 2 || n_k > TB_MAX_GRID) {
    return TB_ERR_INVALID;
  }

  double *flux = malloc((size_t)n_k * (size_t)n_k * sizeof *flux);
  if (!flux) {
    return TB_ERR_MEMORY;
  }

  tb_status_t st =
      tb_berry_flux_grid(model, first_band, n_bands, n_k, flux, NULL);
  if (st == TB_OK) {
    double sum = 0.0;

    for (size_t p = 0; p < (size_t)n_k * (size_t)n_k; p++) {
      sum += flux[p];
    }

    *chern_out = sum / (2.0 * M_PI);
  }

  free(flux);

  return st;
}

// Wilson phase along k1 at fixed k2 using caller-provided scratch:
// blocks holds n_k * n * count entries, work holds count * count
static tb_status_t tb_wilson_core(const tb_model_t *m, int first, int count,
                                  double k2, int n_k, complex_t *blocks,
                                  complex_t *work, double *phase) {
  int n = m->n_orb;
  size_t blk = (size_t)n * (size_t)count;
  double gap_min = INFINITY;
  for (int a = 0; a < n_k; a++) {
    double gap = INFINITY;
    tb_status_t st = tb_group_at(m, 2.0 * M_PI * (double)a / (double)n_k, k2,
                                 first, count, blocks + (size_t)a * blk, &gap);
    if (st != TB_OK) {
      return st;
    }

    gap_min = fmin(gap_min, gap);
  }

  if (!(gap_min >= TB_GAP_TOL)) {
    return TB_ERR_DEGENERATE;
  }

  int bad = 0;
  complex_t prod = c_one();
  for (int a = 0; a < n_k; a++) {
    int a1 = (a + 1) % n_k;

    prod = c_mul(prod, tb_link(n, count, blocks + (size_t)a * blk,
                               blocks + (size_t)a1 * blk, work, &bad));
  }

  if (bad) {
    return TB_ERR_DEGENERATE;
  }

  double g = -atan2(prod.im, prod.re);
  if (g <= -M_PI) {
    g += 2.0 * M_PI;
  }

  *phase = g;

  return TB_OK;
}

tb_status_t tb_wilson_phase(const tb_model_t *model, int first_band,
                            int n_bands, double k_other, int n_k,
                            double *phase_out) {
  if (!tb_group_ok(model, first_band, n_bands) || !phase_out || n_k < 2 ||
      n_k > TB_MAX_GRID || !isfinite(k_other)) {
    return TB_ERR_INVALID;
  }

  size_t blk = (size_t)model->n_orb * (size_t)n_bands;
  if ((size_t)n_k * blk > TB_MAX_STORED_ELEMENTS) {
    return TB_ERR_INVALID;
  }

  complex_t *blocks = malloc((size_t)n_k * blk * sizeof *blocks);
  complex_t *work = malloc((size_t)n_bands * (size_t)n_bands * sizeof *work);
  if (!blocks || !work) {
    free(blocks);
    free(work);

    return TB_ERR_MEMORY;
  }

  tb_status_t st = tb_wilson_core(model, first_band, n_bands, k_other, n_k,
                                  blocks, work, phase_out);

  free(blocks);
  free(work);

  return st;
}

tb_status_t tb_wilson_winding(const tb_model_t *model, int first_band,
                              int n_bands, int n_k1, int n_k2,
                              double *winding_out) {
  if (!tb_group_ok(model, first_band, n_bands) || model->dim != 2 ||
      !winding_out || n_k1 < 2 || n_k1 > TB_MAX_GRID || n_k2 < 4 ||
      n_k2 > TB_MAX_GRID) {
    return TB_ERR_INVALID;
  }

  size_t blk = (size_t)model->n_orb * (size_t)n_bands;
  if ((size_t)n_k1 * blk > TB_MAX_STORED_ELEMENTS) {
    return TB_ERR_INVALID;
  }

  complex_t *blocks = malloc((size_t)n_k1 * blk * sizeof *blocks);
  complex_t *work = malloc((size_t)n_bands * (size_t)n_bands * sizeof *work);
  if (!blocks || !work) {
    free(blocks);
    free(work);

    return TB_ERR_MEMORY;
  }

  tb_status_t st = TB_OK;
  double total = 0.0;
  double prev = 0.0;
  for (int j = 0; j <= n_k2 && st == TB_OK; j++) {
    double cur = 0.0;
    st = tb_wilson_core(model, first_band, n_bands,
                        2.0 * M_PI * (double)j / (double)n_k2, n_k1, blocks,
                        work, &cur);
    if (st != TB_OK) {
      break;
    }

    if (j > 0) {
      double d = cur - prev;

      d -= 2.0 * M_PI * floor((d + M_PI) / (2.0 * M_PI)); // wrap to [-\pi, \pi)
      if (fabs(d) > 0.9 * M_PI) {
        st = TB_ERR_DEGENERATE;

        break;
      }

      total += d;
    }

    prev = cur;
  }

  if (st == TB_OK) {
    *winding_out = total / (2.0 * M_PI);
  }

  free(blocks);
  free(work);

  return st;
}

// Ready-made models
// Add all hoppings from a table; frees model and returns NULL on failure
typedef struct {
  int i;
  int j;
  int n1;
  int n2;
} tb_bond_t;

static tb_model_t *tb_finish(tb_model_t *m, const tb_bond_t *bonds,
                             size_t count, complex_t t) {
  if (!m) {
    return NULL;
  }

  for (size_t b = 0; b < count; b++) {
    if (tb_add_hopping(m, bonds[b].i, bonds[b].j, bonds[b].n1, bonds[b].n2,
                       t) != TB_OK) {
      tb_model_free(m);

      return NULL;
    }
  }

  return m;
}

tb_model_t *tb_model_chain(double t) {
  if (!isfinite(t)) {
    return NULL;
  }

  static const tb_bond_t bonds[] = {{0, 0, 1, 0}};

  return tb_finish(tb_model_alloc(1, 1), bonds, 1, c_real(t));
}

tb_model_t *tb_model_ssh(double v, double w) {
  if (!isfinite(v) || !isfinite(w)) {
    return NULL;
  }

  tb_model_t *m = tb_model_alloc(1, 2);
  if (!m) {
    return NULL;
  }

  if (tb_add_hopping(m, 0, 1, 0, 0, c_real(v)) != TB_OK ||
      tb_add_hopping(m, 1, 0, 1, 0, c_real(w)) != TB_OK) {
    tb_model_free(m);

    return NULL;
  }

  return m;
}

tb_model_t *tb_model_square(double t) {
  if (!isfinite(t)) {
    return NULL;
  }

  static const tb_bond_t bonds[] = {{0, 0, 1, 0}, {0, 0, 0, 1}};

  return tb_finish(tb_model_alloc(2, 1), bonds, 2, c_real(t));
}

tb_model_t *tb_model_graphene(double t) {
  if (!isfinite(t)) {
    return NULL;
  }

  static const tb_bond_t bonds[] = {{0, 1, 0, 0}, {0, 1, -1, 0}, {0, 1, 0, -1}};

  return tb_finish(tb_model_alloc(2, 2), bonds, 3, c_real(t));
}

tb_model_t *tb_model_haldane(double t1, double t2, double phi, double mass) {
  if (!isfinite(t1) || !isfinite(t2) || !isfinite(phi) || !isfinite(mass)) {
    return NULL;
  }

  tb_model_t *m = tb_model_graphene(t1);
  if (!m) {
    return NULL;
  }

  static const int nu[3][2] = {{1, 0}, {-1, 1}, {0, -1}};
  complex_t ta = c_from_polar(t2, phi);
  complex_t tb = c_from_polar(t2, -phi);
  int ok =
      tb_set_onsite(m, 0, mass) == TB_OK && tb_set_onsite(m, 1, -mass) == TB_OK;
  for (int k = 0; k < 3 && ok; k++) {
    ok = tb_add_hopping(m, 0, 0, nu[k][0], nu[k][1], ta) == TB_OK &&
         tb_add_hopping(m, 1, 1, nu[k][0], nu[k][1], tb) == TB_OK;
  }

  if (!ok) {
    tb_model_free(m);

    return NULL;
  }

  return m;
}

tb_model_t *tb_model_kagome(double t) {
  if (!isfinite(t)) {
    return NULL;
  }

  static const tb_bond_t bonds[] = {{0, 1, 0, 0}, {0, 1, -1, 0},
                                    {0, 2, 0, 0}, {0, 2, 0, -1},
                                    {1, 2, 0, 0}, {1, 2, 1, -1}};

  return tb_finish(tb_model_alloc(2, 3), bonds, 6, c_real(t));
}

tb_model_t *tb_model_hofstadter(int p, int q, double t) {
  if (p < 1 || q < 2 || p >= q || q > TB_MAX_ORBITALS || !isfinite(t)) {
    return NULL;
  }

  tb_model_t *m = tb_model_alloc(2, q);
  if (!m) {
    return NULL;
  }

  double alpha = (double)p / (double)q;
  int ok = 1;
  for (int i = 0; i < q - 1 && ok; i++) {
    ok = tb_add_hopping(m, i, i + 1, 0, 0, c_real(-t)) == TB_OK;
  }

  ok = ok && tb_add_hopping(m, q - 1, 0, -1, 0, c_real(-t)) == TB_OK;
  for (int i = 0; i < q && ok; i++) {
    ok = tb_add_hopping(m, i, i, 0, 1,
                        c_from_polar(t, 2.0 * M_PI * alpha * (double)i)) ==
         TB_OK;
  }

  if (!ok) {
    tb_model_free(m);

    return NULL;
  }

  return m;
}

// Analytic references
double tb_chain_energy(double t, double k) { return 2.0 * t * cos(k); }

double tb_square_energy(double t, double k1, double k2) {
  return 2.0 * t * (cos(k1) + cos(k2));
}

double tb_graphene_energy(double t, double k1, double k2) {
  double re = 1.0 + cos(k1) + cos(k2);
  double im = sin(k1) + sin(k2);

  return fabs(t) * hypot(re, im);
}

double tb_haldane_valley_gap(double t2, double phi, double mass, int valley) {
  if ((valley != 1 && valley != -1) || !isfinite(t2) || !isfinite(phi) ||
      !isfinite(mass)) {
    return NAN;
  }

  return 2.0 * fabs(mass - (double)valley * 3.0 * sqrt(3.0) * t2 * sin(phi));
}

int tb_haldane_chern(double t1, double t2, double phi, double mass) {
  if (!isfinite(t1) || !isfinite(t2) || !isfinite(phi) || !isfinite(mass) ||
      t1 == 0.0 || fabs(t2) > fabs(t1) / 3.0) {
    return TB_CHERN_UNDEFINED;
  }

  double drive = t2 * sin(phi);
  double limit = 3.0 * sqrt(3.0) * fabs(drive);
  if (fabs(fabs(mass) - limit) <= 1e-12 * fmax(1.0, limit)) {
    return TB_CHERN_UNDEFINED;
  }

  if (fabs(mass) < limit) {
    return drive > 0.0 ? 1 : -1;
  }

  return 0;
}

double tb_ssh_zak_phase(double v, double w) {
  if (!isfinite(v) || !isfinite(w) || fabs(v) == fabs(w)) {
    return NAN;
  }

  return fabs(v) > fabs(w) ? 0.0 : M_PI;
}

double tb_kagome_flat_band_energy(double t) { return -2.0 * t; }
