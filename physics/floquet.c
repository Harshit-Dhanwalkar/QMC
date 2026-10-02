#include "floquet.h"
#include "../core/complex.h"
#include "../core/linalg/complex_eigh.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Dimension limits keep n^3 / (n (2K+1))^2 work and memory sane and prevent
// integer overflow in size computations
#define FLOQUET_MAX_DIM 512
#define FLOQUET_MAX_SAMBE_DIM 1500

static double floquet_fold(double eps, double omega) {
  double x = fmod(eps + 0.5 * omega, omega);
  if (x < 0.0) {
    x += omega;
  }

  return x - 0.5 * omega;
}

static floquet_result_t *floquet_result_alloc(int n, double omega) {
  floquet_result_t *res = calloc(1, sizeof *res);
  if (!res) {
    return NULL;
  }

  res->n = n;
  res->omega = omega;
  res->quasienergies = calloc((size_t)n, sizeof *res->quasienergies);
  res->modes = cmatrix_alloc(n, n);
  if (!res->quasienergies || !res->modes) {
    floquet_result_free(res);

    return NULL;
  }

  return res;
}

void floquet_result_free(floquet_result_t *result) {
  if (!result) {
    return;
  }

  free(result->quasienergies);
  cmatrix_free(result->modes);
  free(result);
}

/* Time-domain solver */
// out = -i * H * X for n x n row-major arrays (out must not alias X)
static void floquet_neg_i_hx(int n, const complex_t *H, const complex_t *X,
                             complex_t *out) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double re = 0.0;
      double im = 0.0;
      for (int k = 0; k < n; k++) {
        complex_t h = H[(size_t)i * n + k];
        complex_t x = X[(size_t)k * n + j];
        re += h.re * x.re - h.im * x.im;
        im += h.re * x.im + h.im * x.re;
      }

      // multiply by -i: (re + i im) * (-i) = im - i re
      out[(size_t)i * n + j] = (complex_t){im, -re};
    }
  }
}

// dst = U + scale * k
static void floquet_axpy(size_t len, const complex_t *U, double scale,
                         const complex_t *k, complex_t *dst) {
  for (size_t i = 0; i < len; i++) {
    dst[i].re = U[i].re + scale * k[i].re;
    dst[i].im = U[i].im + scale * k[i].im;
  }
}

/* Propagate U(T) with RK4 on dU/dt = -i H(t) U, U(0) = I. U is n x n row-major
 *
 * Returns 0 on success, -1 on allocation failure
 */
static int floquet_propagate(int n, floquet_hamiltonian_fn Hfn, void *params,
                             double T, int steps, complex_t *U) {
  size_t len = (size_t)n * (size_t)n;
  complex_t *work = calloc(6 * len, sizeof *work);
  cmatrix_t *Ha = cmatrix_alloc(n, n);
  cmatrix_t *Hm = cmatrix_alloc(n, n);
  cmatrix_t *Hb = cmatrix_alloc(n, n);
  if (!work || !Ha || !Hm || !Hb) {
    free(work);
    cmatrix_free(Ha);
    cmatrix_free(Hm);
    cmatrix_free(Hb);

    return -1;
  }

  complex_t *k1 = work;
  complex_t *k2 = work + len;
  complex_t *k3 = work + 2 * len;
  complex_t *k4 = work + 3 * len;
  complex_t *tmp = work + 4 * len;

  for (size_t i = 0; i < len; i++) {
    U[i] = c_zero();
  }
  for (int i = 0; i < n; i++) {
    U[(size_t)i * n + i] = c_one();
  }

  double h = T / (double)steps;
  memset(Ha->data, 0, len * sizeof *Ha->data);
  Hfn(0.0, params, Ha);
  for (int s = 0; s < steps; s++) {
    double t = (double)s * h;
    memset(Hm->data, 0, len * sizeof *Hm->data);
    memset(Hb->data, 0, len * sizeof *Hb->data);

    Hfn(t + 0.5 * h, params, Hm);
    Hfn((double)(s + 1) * h, params, Hb);

    floquet_neg_i_hx(n, Ha->data, U, k1);
    floquet_axpy(len, U, 0.5 * h, k1, tmp);
    floquet_neg_i_hx(n, Hm->data, tmp, k2);
    floquet_axpy(len, U, 0.5 * h, k2, tmp);
    floquet_neg_i_hx(n, Hm->data, tmp, k3);
    floquet_axpy(len, U, h, k3, tmp);
    floquet_neg_i_hx(n, Hb->data, tmp, k4);

    for (size_t i = 0; i < len; i++) {
      U[i].re +=
          h / 6.0 * (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re);
      U[i].im +=
          h / 6.0 * (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im);
    }

    // H at end of this step is H at start of next one
    cmatrix_t *swap = Ha;
    Ha = Hb;
    Hb = swap;
  }

  free(work);
  cmatrix_free(Ha);
  cmatrix_free(Hm);
  cmatrix_free(Hb);

  return 0;
}

static int floquet_all_finite(const complex_t *a, size_t len) {
  for (size_t i = 0; i < len; i++) {
    if (!isfinite(a[i].re) || !isfinite(a[i].im)) {
      return 0;
    }
  }

  return 1;
}

// max |U^\dagger U - I|
static double floquet_unitarity_error(int n, const complex_t *U) {
  double err = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double re = 0.0;
      double im = 0.0;

      for (int k = 0; k < n; k++) {
        complex_t a = U[(size_t)k * n + i]; // (U^\dagger)_{ik} = conj(U_{ki})
        complex_t b = U[(size_t)k * n + j];

        re += a.re * b.re + a.im * b.im;
        im += a.re * b.im - a.im * b.re;
      }

      if (i == j) {
        re -= 1.0;
      }

      err = fmax(err, hypot(re, im));
    }
  }

  return err;
}

/* Diagonalize unitary U via Hermitian A = C + s S
 *
 * Where
 *    C = (U + U^\dagger)/2 and S = (U - U^\dagger)/(2i) commute
 *
 * Fills `\lambda` (n eigenvalues of U) and `vecs` (n x n, columns =
 * eigenvectors) from attempt with lowest residual; returns that residual (or
 * \infty if every attempt failed outright)
 */
static double floquet_unitary_eig(int n, const complex_t *U, complex_t *lambda,
                                  cmatrix_t *vecs) {
  static const double shifts[] = {0.6180339887498949, 1.3247179572447460,
                                  -0.7548776662466927, 2.2360679774997897,
                                  -1.4142135623730951};
  double best = INFINITY;
  cmatrix_t *A = cmatrix_alloc(n, n);
  if (!A) {
    return best;
  }

  complex_t *lam = malloc((size_t)n * sizeof *lam);
  if (!lam) {
    cmatrix_free(A);

    return best;
  }

  for (size_t attempt = 0; attempt < sizeof shifts / sizeof shifts[0];
       attempt++) {
    double s = shifts[attempt];

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        complex_t uij = U[(size_t)i * n + j];
        complex_t uji_c = c_conj(U[(size_t)j * n + i]);
        complex_t sum = c_add(uij, uji_c); // 2 C_ij
        complex_t dif = c_sub(uij, uji_c); // 2i S_ij
        // A_ij = \sum/2 + s * dif/(2i) = \sum/2 - i s dif / 2
        complex_t term = c_new(0.5 * sum.re + 0.5 * s * dif.im,
                               0.5 * sum.im - 0.5 * s * dif.re);

        A->data[(size_t)i * n + j] = term;
      }
    }

    eigen_t *eig = cmatrix_eigh_complex(A);
    if (!eig) {
      continue;
    }

    double residual = 0.0;
    for (int k = 0; k < n; k++) {
      // \lambda_k = <v|U|v>
      complex_t rq = c_zero();
      for (int i = 0; i < n; i++) {
        complex_t uv = c_zero();
        for (int j = 0; j < n; j++) {
          uv = c_add(uv, c_mul(U[(size_t)i * n + j],
                               eig->eigenvectors->data[(size_t)j * n + k]));
        }

        rq = c_add(
            rq, c_mul(c_conj(eig->eigenvectors->data[(size_t)i * n + k]), uv));
      }

      lam[k] = rq;
      double r2 = 0.0;
      for (int i = 0; i < n; i++) {
        complex_t uv = c_zero();

        for (int j = 0; j < n; j++) {
          uv = c_add(uv, c_mul(U[(size_t)i * n + j],
                               eig->eigenvectors->data[(size_t)j * n + k]));
        }

        complex_t d =
            c_sub(uv, c_mul(rq, eig->eigenvectors->data[(size_t)i * n + k]));
        r2 += c_abs2(d);
      }

      residual = fmax(residual, sqrt(r2));
    }

    if (residual < best) {
      best = residual;

      memcpy(lambda, lam, (size_t)n * sizeof *lam);
      memcpy(vecs->data, eig->eigenvectors->data,
             (size_t)n * (size_t)n * sizeof *vecs->data);
    }

    eigen_free(eig);

    if (best < 1e-8) {
      break;
    }
  }

  free(lam);
  cmatrix_free(A);

  return best;
}

// Sort quasi-energies ascending, permuting columns of `modes` alongside.
static void floquet_sort(int n, double *eps, cmatrix_t *modes) {
  for (int i = 1; i < n; i++) {
    double e = eps[i];
    int j = i - 1;
    // insertion sort with column swaps
    int pos = i;
    while (j >= 0 && eps[j] > e) {
      eps[j + 1] = eps[j];

      for (int r = 0; r < n; r++) {
        complex_t tmp = modes->data[(size_t)r * n + j + 1];

        modes->data[(size_t)r * n + j + 1] = modes->data[(size_t)r * n + j];
        modes->data[(size_t)r * n + j] = tmp;
      }

      pos = j;
      j--;
    }

    eps[pos] = e;
  }
}

floquet_result_t *floquet_solve_time(int n, floquet_hamiltonian_fn H,
                                     void *params, double omega, int steps) {
  if (n < 1 || n > FLOQUET_MAX_DIM || !H || !(omega > 0.0) ||
      !isfinite(omega) || steps < 4) {
    return NULL;
  }

  size_t len = (size_t)n * (size_t)n;
  complex_t *U = malloc(len * sizeof *U);
  complex_t *lambda = calloc((size_t)n, sizeof *lambda);
  floquet_result_t *res = floquet_result_alloc(n, omega);
  if (!U || !lambda || !res) {
    free(U);
    free(lambda);
    floquet_result_free(res);

    return NULL;
  }

  double T = 2.0 * M_PI / omega;
  if (floquet_propagate(n, H, params, T, steps, U) != 0 ||
      !floquet_all_finite(U, len)) {
    free(U);
    free(lambda);
    floquet_result_free(res);

    return NULL;
  }

  res->unitarity_error = floquet_unitarity_error(n, U);
  res->eigen_residual = floquet_unitary_eig(n, U, lambda, res->modes);
  if (!isfinite(res->eigen_residual)) {
    free(U);
    free(lambda);
    floquet_result_free(res);

    return NULL;
  }

  for (int k = 0; k < n; k++) {
    double phase = atan2(lambda[k].im, lambda[k].re); // (-\pi, \pi]
    double eps = -phase / T;
    if (eps >= 0.5 * omega) {
      eps -= omega;
    }

    res->quasienergies[k] = eps;
  }

  floquet_sort(n, res->quasienergies, res->modes);

  free(U);
  free(lambda);

  return res;
}

/* Sambe-space solver */
floquet_result_t *floquet_solve_sambe(int n, cmatrix_t *const *harmonics, int M,
                                      double omega, int K) {
  if (n < 1 || n > FLOQUET_MAX_DIM || !harmonics || M < 0 || K < M + 1 ||
      !(omega > 0.0) || !isfinite(omega)) {
    return NULL;
  }

  long dim_l = (long)n * (2L * K + 1L);
  if (dim_l > FLOQUET_MAX_SAMBE_DIM) {
    return NULL;
  }

  int dim = (int)dim_l;

  // Validate harmonics: shapes, finiteness, and H_{-m} = H_m^\dagger
  double scale = 1.0;
  for (int m = -M; m <= M; m++) {
    const cmatrix_t *h = harmonics[m + M];
    if (!h || h->nrows != n || h->ncols != n ||
        !floquet_all_finite(h->data, (size_t)n * (size_t)n)) {
      return NULL;
    }

    for (size_t i = 0; i < (size_t)n * (size_t)n; i++) {
      scale = fmax(scale, hypot(h->data[i].re, h->data[i].im));
    }
  }

  for (int m = 0; m <= M; m++) {
    const cmatrix_t *hp = harmonics[M + m];
    const cmatrix_t *hm = harmonics[M - m];
    for (int a = 0; a < n; a++) {
      for (int b = 0; b < n; b++) {
        complex_t d = c_sub(hm->data[(size_t)a * n + b],
                            c_conj(hp->data[(size_t)b * n + a]));
        if (hypot(d.re, d.im) > 1e-9 * scale) {
          return NULL;
        }
      }
    }
  }

  cmatrix_t *HF = cmatrix_alloc(dim, dim);
  floquet_result_t *res = floquet_result_alloc(n, omega);
  if (!HF || !res) {
    cmatrix_free(HF);
    floquet_result_free(res);

    return NULL;
  }

  for (int a = -K; a <= K; a++) {
    for (int b = -K; b <= K; b++) {
      int m = a - b;
      if (m < -M || m > M) {
        continue;
      }

      const cmatrix_t *h = harmonics[m + M];
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          HF->data[(size_t)((a + K) * n + i) * dim +
                   (size_t)((b + K) * n + j)] = h->data[(size_t)i * n + j];
        }
      }
    }

    for (int i = 0; i < n; i++) {
      size_t idx = (size_t)((a + K) * n + i) * dim + (size_t)((a + K) * n + i);
      HF->data[idx].re += (double)a * omega;
    }
  }

  eigen_t *eig = cmatrix_eigh_complex(HF);

  cmatrix_free(HF);
  if (!eig) {
    floquet_result_free(res);

    return NULL;
  }

  // Central n eigenvalues: one replica of every quasi-energy class.
  int first = K * n;
  int n_blocks = 2 * K + 1;
  double worst_edge = 0.0;
  for (int j = 0; j < n; j++) {
    int col = first + j;
    res->quasienergies[j] = floquet_fold(eig->eigenvalues[col], omega);

    double edge = 0.0;
    for (int blk = 0; blk < n_blocks; blk++) {
      for (int i = 0; i < n; i++) {
        complex_t v =
            eig->eigenvectors->data[(size_t)(blk * n + i) * dim + (size_t)col];
        res->modes->data[(size_t)i * n + j] =
            c_add(res->modes->data[(size_t)i * n + j], v);
        if (blk == 0 || blk == n_blocks - 1) {
          edge += c_abs2(v);
        }
      }
    }

    worst_edge = fmax(worst_edge, edge);
  }

  res->truncation_weight = worst_edge;

  eigen_free(eig);

  floquet_sort(n, res->quasienergies, res->modes);

  return res;
}

/* Stroboscopic evolution and analytic references */
cvector_t *floquet_stroboscopic_state(const floquet_result_t *result,
                                      const cvector_t *psi0, long n_periods) {
  if (!result || !psi0 || !psi0->data || psi0->n != result->n ||
      n_periods < 0) {
    return NULL;
  }

  int n = result->n;
  cvector_t *out = cvector_alloc(n);
  if (!out) {
    return NULL;
  }

  double T = 2.0 * M_PI / result->omega;
  for (int i = 0; i < n; i++) {
    out->data[i] = c_zero();
  }

  for (int a = 0; a < n; a++) {
    complex_t c = c_zero(); // <u_a | \psi0>
    for (int i = 0; i < n; i++) {
      c = c_add(c, c_mul(c_conj(result->modes->data[(size_t)i * n + a]),
                         psi0->data[i]));
    }

    double angle =
        fmod(-result->quasienergies[a] * T * (double)n_periods, 2.0 * M_PI);
    complex_t w = c_mul(c, c_from_polar(1.0, angle));
    for (int i = 0; i < n; i++) {
      out->data[i] =
          c_add(out->data[i], c_mul(w, result->modes->data[(size_t)i * n + a]));
    }
  }

  return out;
}

double floquet_circular_rabi_quasienergy(double omega0, double Omega,
                                         double omega, int sign) {
  if (!isfinite(omega0) || !isfinite(Omega) || !(omega > 0.0) ||
      !isfinite(omega) || (sign != 1 && sign != -1)) {
    return NAN;
  }

  double d = omega0 - omega;
  double omega_r = sqrt(d * d + Omega * Omega);

  return floquet_fold(0.5 * omega + 0.5 * (double)sign * omega_r, omega);
}

double floquet_bessel_j0(double x) {
  if (!isfinite(x) || fabs(x) > 1000.0) {
    return NAN;
  }

  // NOTE: J0(x) = (1/\pi) \int_0^\pi \cos(x \sin t) dt; integrand is smooth and
  // periodic, so trapezoid rule converges geometrically once point count
  // exceeds ~|x|
  int pts = 64 + 2 * (int)fabs(x);
  double h = M_PI / (double)pts;
  double sum = 0.5 * (cos(0.0) + cos(x * sin(M_PI)));
  for (int k = 1; k < pts; k++) {
    sum += cos(x * sin((double)k * h));
  }

  return sum * h / M_PI;
}

double floquet_cdt_splitting(double Delta, double A, double omega) {
  if (!isfinite(Delta) || !isfinite(A) || !(omega > 0.0) || !isfinite(omega)) {
    return NAN;
  }

  return fabs(Delta) * fabs(floquet_bessel_j0(A / omega));
}
