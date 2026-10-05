/*
 * qmc_butterfly.c - WebAssembly front end for the Hofstadter butterfly.
 *
 * Everything physical comes from the library's tight-binding module
 * (physics/tight_binding.c): tb_model_hofstadter() builds the q-orbital
 * magnetic unit cell, tb_bands() diagonalises it, tb_chern_number() computes
 * the Fukui-Hatsugai-Suzuki Chern number of the filled bands. This file only
 * packages three questions for JavaScript:
 *
 *   qmc_bf_edges   - the q band intervals [lo, hi] at flux p/q
 *   qmc_bf_tknn    - Hall number of gap r from the TKNN Diophantine equation
 *   qmc_bf_chern   - the same number computed numerically by the library
 *
 * Band edges: for H(k1,k2) of tb_model_hofstadter the characteristic
 * polynomial depends on k1 and k2 only through cos(k1) and cos(q*k2), so
 * every band edge sits at one of k1 in {0, pi}, k2 in {0, pi/q}. Four
 * diagonalisations per flux therefore give the exact band intervals (checked
 * against a dense k-grid for all coprime p/q with q <= 14 in
 * wasm/test/butterfly_check.c).
 */

#include "../physics/tight_binding.h"
#include <limits.h>
#include <math.h>
#include <stddef.h>

#ifdef __wasm__
#define QMC_EXPORT(name) __attribute__((export_name(#name)))
#else
#define QMC_EXPORT(name) /* native builds (tests) need no export attribute */
#endif

enum { BF_MAX_Q = TB_MAX_ORBITALS };

static double g_buf[2 * BF_MAX_Q]; /* scratch shared with JavaScript */

/* Pointer to the scratch buffer (2*q doubles: lo0, hi0, lo1, hi1, ...). */
QMC_EXPORT(qmc_bf_buffer)
double *qmc_bf_buffer(void) { return g_buf; }

QMC_EXPORT(qmc_bf_max_q)
int qmc_bf_max_q(void) { return BF_MAX_Q; }

static int valid_flux(int p, int q) {
  return q >= 2 && q <= BF_MAX_Q && p >= 1 && p < q;
}

/* Band intervals at flux p/q (hopping t = 1) written to the scratch buffer as
 * (lo, hi) pairs in ascending band order. Returns 0 on success. */
QMC_EXPORT(qmc_bf_edges)
int qmc_bf_edges(int p, int q) {
  if (!valid_flux(p, q)) {
    return -1;
  }
  tb_model_t *m = tb_model_hofstadter(p, q, 1.0);
  if (!m) {
    return -2;
  }
  for (int b = 0; b < q; b++) {
    g_buf[2 * b] = INFINITY;
    g_buf[2 * b + 1] = -INFINITY;
  }
  const double k1s[2] = {0.0, M_PI};
  const double k2s[2] = {0.0, M_PI / q};
  double e[BF_MAX_Q];
  int rc = 0;
  for (int i = 0; i < 2 && rc == 0; i++) {
    for (int j = 0; j < 2; j++) {
      if (tb_bands(m, k1s[i], k2s[j], e) != TB_OK) {
        rc = -3;
        break;
      }
      for (int b = 0; b < q; b++) {
        if (e[b] < g_buf[2 * b]) {
          g_buf[2 * b] = e[b];
        }
        if (e[b] > g_buf[2 * b + 1]) {
          g_buf[2 * b + 1] = e[b];
        }
      }
    }
  }
  tb_model_free(m);
  return rc;
}

static int inv_mod(int a, int m) {
  for (int x = 1; x < m; x++) {
    if ((a * x) % m == 1) {
      return x;
    }
  }
  return 0;
}

/* Hall number t_r of the gap above the lowest r bands (r = 1 .. q-1), from the
 * TKNN Diophantine equation r = q*s + p*t with |t| <= q/2. For even q the
 * middle gap (r = q/2) closes at the Dirac points and has no Hall number:
 * returns INT_MIN there, and for invalid input or non-coprime p, q. */
QMC_EXPORT(qmc_bf_tknn)
int qmc_bf_tknn(int p, int q, int r) {
  if (!valid_flux(p, q) || r < 1 || r >= q) {
    return INT_MIN;
  }
  int inv = inv_mod(p % q, q);
  if (inv == 0) {
    return INT_MIN; /* p, q not coprime */
  }
  if (q % 2 == 0 && 2 * r == q) {
    return INT_MIN;
  }
  int t = (r * inv) % q;
  if (2 * t > q) {
    t -= q;
  }
  return t;
}

/* Chern number of the lowest r bands at flux p/q, computed numerically by the
 * library on an n_k x n_k grid (Fukui-Hatsugai-Suzuki). Returns NaN when the
 * gap is closed or the grid is too coarse. Also usable to cross-check
 * qmc_bf_tknn: the two agree for every gapped r. */
QMC_EXPORT(qmc_bf_chern)
double qmc_bf_chern(int p, int q, int r, int n_k) {
  if (!valid_flux(p, q) || r < 1 || r >= q || n_k < 4) {
    return NAN;
  }
  tb_model_t *m = tb_model_hofstadter(p, q, 1.0);
  if (!m) {
    return NAN;
  }
  double c = NAN;
  if (tb_chern_number(m, 0, r, n_k, &c) != TB_OK) {
    c = NAN;
  }
  tb_model_free(m);
  return c;
}
